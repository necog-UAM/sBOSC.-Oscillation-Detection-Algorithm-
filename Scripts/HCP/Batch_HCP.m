%% Batch script to run sBOSC in HCP data. 

% sBOSC path
p.sBOSC = 'Z:\Toolbox\sourceBOSC';
addpath(genpath(p.sBOSC))

% Fieldtrip path
addpath('Z:\Toolbox\fieldtrip-20230118')  % Fieldtrip path
ft_defaults

% OMEGA data path
p.data = 'Z:\HCP\data\'; % data results
p.raw = 'Z:\HCP\raw_preproc\'; % preprocessed data      
p.scripts = 'Z:\OMEGA\Enrique\OMEGA-Repo\Scripts\HCP';

cd (p.data)
temp = ls;
subs = deblank(temp(3:end,:));                      % all but  .  and  ..  
tasks    = {'Motort'};
sessions = [10 11];
Ncondit  = 3;

Nsub = size(subs,1);
substsk = {};
task = tasks{1};
ct = 1;
for s = 1: Nsub
    sub = deblank(subs(s,:));
    cd([p.raw sub '\MEG'])
    if exist(tasks{1})==7
        substsk{1}{ct} = sub;
        if s==58 
            ct = ct;             % to remove subject 36 in Motort (only 30 trials in one of the conditions)
        else
            ct = ct+1;
        end
    end
end

% Motort
condit_label{1} = 'LeftHand';
condit_label{2} = 'RightHand';

%% ----------- 1. Reading, preprocessing and artifact correction ---------------- %
% The HCP dataset has already been preprocessed

%% sourceBOSC
%% ----------------------- 2. Source reconstruction  ---------------------------- %

% MEG-MRI coregistration: mark fiducials (lpa, rpa and nasion; then quit) 
% and check result (if coregistration is not accurate, repeat procedure)
for s = 1:Nsub
    sub = deblank(subs(s,:));
    sprintf(['\n Processing Sub:' sub '..\n'])
    task = [tasks{1} '11'];
    cd([p.raw sub '\MEG'])
    datas = {};
    if exist(tasks{1})==7
        mri = ft_read_mri([p.raw sub '\MEG\anatomy\T1w_acpc_dc_restore.nii']);        
        fid  = fopen([p.raw sub '\MEG\anatomy\' num2str(sub) '_MEG_anatomy_transform.txt']);        % coregistration: load transformation matrix from HCP database
        txt  = textscan(fid,'%s');
        txt  = txt{1};
        fclose(fid);

        for tt = 1:length(txt)
            if strcmp(txt{tt},'transform.vox07mm2bti')
                t1 = tt;
            end
        end

        transform = [];
        for i = 1:16
            transform(i) = str2num(txt{t1+i+2});
        end
        mri.transform = reshape(transform,[4 4])';

        for ses = sessions(1):sessions(end)
            cd([p.raw sub '\MEG\Motort\tmegpreproc'])
            a = ls('*Motort_tmegpreproc_TEMG*');
            for i = 1:size(a,1)
                temp = strfind(a(i,:),['MEG_' num2str(ses) '-Motort_tmegpreproc_TEMG']);
                if temp > 1
                    disp(['Reading ' a(i,:)])
                    load(a(i,:))
                end
            end
 
            % Remove EMG channels in motor task
            cfg         = [];
            cfg.channel = 'MEG';              
            data        = ft_preprocessing(cfg,data);

            % Remove NaN trials
            trls = [];
            trlnan = [];
            for i = 1:length(data.trial)
                trls(i) = size(data.trial{i},2);
                trlnan(i) = sum(sum(isnan(data.trial{i})));
            end

            % Remove shorter trials
            cfg = [];
            cfg.trials = trls>=610 & trlnan==0;           % remove shorter trials and trials with NaNs
            data       = ft_redefinetrial(cfg,data);
    
            % Select condition
            ct = 1;
            for cond = [1 4]            % 1-Left Hand,  2 - Left Foot, 4 - Right Hand. 5 - Right Foot, 6 - Fixation
                cfg = [];
                cfg.trials = data.trialinfo(:,2) == cond;    
                datacond   = ft_redefinetrial(cfg,data);
    
                cfg = [];
                cfg.mri = mri;
                cfg.figure = 'yes';     % to check that sensors, volyume and grid are correctly coregistered
                [source_forward, source_inverse, datasource] = sBOSC_beamforming(datacond, cfg);   % do you want to change the anatomical labels for the axes [Y, n]? Y (r,a,s,i)
                
                datas{ct}{ses} = datasource;
                ct = ct + 1;
            end
        end      % session end

        for cond = 1:length(condit_label)
            % Combine both sessions
            cfg=[];
            data_combined = ft_appenddata(cfg, datas{cond}{sessions(1)}, datas{cond}{sessions(2)});     
        
            % Downsample
            cfg = [];
            cfg.resamplefs = 256;
            data_combined = ft_resampledata(cfg, data_combined);

            % Single to save space
            datasource = ft_struct2single(data_combined);

            % Save
            save([p.data sub '\'  task '\datasource_' condit_label{cond} '.mat'], 'datasource');
        end
    end
end

%% ----------- 3. Compute aperiodic to get thresholds ---------------- %

for s = 1:Nsub
    for cond = 1:length(condit_label)
        sub = deblank(subs(s,:));
        task = [tasks{1} '11'];
        cd([p.data '\' sub])
        if exist(task)==7
            sprintf(['\n Processing Sub: ' num2str(s) '/' num2str(Nsub) '..\n'])
            load([p.data sub '\'  task '\datasource_' condit_label{cond} '.mat'])
           
                cfg = [];
                cfg.datatype   = 'trials';
            evalc('aperiodic = sBOSC_Aperiodic(datasource, cfg)');

            % Save
            save([p.data sub  '\' task '\aperiodic_' condit_label{cond} '.mat'], aperiodic)
        end
    end
end

%% ----------- 4. Power-spectrum and power threshold ---------------- %
for s = 1:Nsub 
    sprintf(['\n Processing Sub:' subs{s} '..\n'])
    load([p.data sub '\'  task '\datasource_' condit_label{cond} '.mat'])
    load([p.data sub  '\' task '\aperiodic_' condit_label{cond} '.mat'])
    
        cfg = [];
        cfg.frex = exp(0.6:0.1:3.7);
        cfg.apthshld = 95;
    evalc('[powspctm, thshld, frex, fsample] = sBOSC_timefreq(datasource, aperiodic, cfg)');

%% ----------- 5. Compute local and spatial peaks ---------------- %
        cfg = [];
    evalc('[spatialpks, localpks] = sBOSC_spatialpeaks(powspctm, thshld, cfg)');

%% ----------- 6. Construct episodes ---------------- %
        cfg = [];
        cfg.frex = frex;
        cfg.fsample = fsample;
        cfg.min_cycles = 1;
    evalc('[episodes, episocc] = sBOSC_episodes(spatialpks, powspctm, cfg)');

%% ----------- 7. Connect episodes ---------------- %

        cfg = [];
        cfg.fsample = fsample;
        cfg.frex = frex;
        cfg.time = size(spatialpks,4);
    evalc('[conepis, conepisocc] = sBOSC_connect_episodes(cfg, episodes)');

end

%% Evaluate
frex = exp(0.6:0.1:3.7); % 1.8 Hz to 40 Hz
nFrex = length(frex);

load(['Z:\HCP\data\140117\Motort11/datasource_RightHand.mat'])
nTp = ceil(length(datasource.time{1})/2);
time = downsample(datasource.time{1}, 2);
t0 =   dsearchn(time',0);

nVox = 1925; 

% Initialize
trial_durcyc_pre  = cell(Nsub, 1);
trial_pow_pre     = cell(Nsub, 1);
trial_durcyc_post = cell(Nsub, 1);
trial_pow_post    = cell(Nsub, 1);
trials_per_sub    = zeros(Nsub, 1);

it = 0;
%% Loop
for s = 1:Nsub
    s
    sub = deblank(subs(s,:));
    task = [tasks{1} '11'];
    sub_task_path = fullfile(p.data, sub, task);
    
    if exist(sub_task_path, 'dir') == 7
        it = it + 1;
        epis_path = fullfile(sub_task_path, 'conepis', 'RightHand');
        
        % Load one vox to check epis size
        load(fullfile(epis_path, 'epis_v1.mat'), 'epis');
        nTrials = length(epis);
        
        % Store number of trials
        trials_per_sub(it) = nTrials; 
        
        % Initialize
        sub_trial_durcyc_pre  = zeros(nTrials, nVox, nFrex, 'single');
        sub_trial_pow_pre     = zeros(nTrials, nVox, nFrex, 'single');
        sub_trial_count_pre   = zeros(nTrials, nVox, nFrex, 'single');
        
        sub_trial_durcyc_post = zeros(nTrials, nVox, nFrex, 'single');
        sub_trial_pow_post    = zeros(nTrials, nVox, nFrex, 'single');
        sub_trial_count_post  = zeros(nTrials, nVox, nFrex, 'single');

        subepis = cell(nTrials, nVox);
        for v = 1:nVox
            load(fullfile(epis_path, ['epis_v' num2str(v) '.mat']), 'epis');
            subepis(:,v) = epis;
        end
        
        for trl = 1:nTrials
            for v = 1:nVox
                tempepis = subepis{trl,v};
                for ep = 1:length(tempepis)
                    
                    % PRE
                    if tempepis(ep).timeps(1) < t0 
                        ftemp = dsearchn(frex', tempepis(ep).freq);
                                       
                        sub_trial_durcyc_pre(trl,v,ftemp) = sub_trial_durcyc_pre(trl,v,ftemp) + tempepis(ep).dur_cyc;
                        sub_trial_pow_pre(trl,v,ftemp)    = sub_trial_pow_pre(trl,v,ftemp) + tempepis(ep).power;
                        sub_trial_count_pre(trl,v,ftemp)  = sub_trial_count_pre(trl,v,ftemp) + 1;
                  
                    % POST
                    elseif tempepis(ep).timeps(1) > t0 
                        ftemp2 = dsearchn(frex', tempepis(ep).freq);
                        
                        sub_trial_durcyc_post(trl,v,ftemp2) = sub_trial_durcyc_post(trl,v,ftemp2) + tempepis(ep).dur_cyc;
                        sub_trial_pow_post(trl,v,ftemp2)    = sub_trial_pow_post(trl,v,ftemp2) + tempepis(ep).power;
                        sub_trial_count_post(trl,v,ftemp2)  = sub_trial_count_post(trl,v,ftemp2) + 1;
                    end
                end
            end
        end
        % Normalization by number of trials (change 0 to 1)
        count_pre  = max(sub_trial_count_pre, 1); 
        count_post = max(sub_trial_count_post, 1);
        
        % Mean power by trial
        mean_trial_pow_pre  = sub_trial_pow_pre  ./ count_pre;
        mean_trial_pow_post = sub_trial_pow_post ./ count_post;

        trial_durcyc_pre{it}  = sub_trial_durcyc_pre;
        trial_pow_pre{it}     = mean_trial_pow_pre;
        
        trial_durcyc_post{it} = sub_trial_durcyc_post;
        trial_pow_post{it}    = mean_trial_pow_post;

    end
end
% Save

% Clear
trial_durcyc_pre  = trial_durcyc_pre(1:it);
trial_pow_pre     = trial_pow_pre(1:it);
trial_durcyc_post = trial_durcyc_post(1:it);
trial_pow_post    = trial_pow_post(1:it);
trials_per_sub    = trials_per_sub(1:it);

save('allsubs_pow-dur_epis_RightHand.mat', 'trial_durcyc_pre', 'trial_pow_pre', 'trial_durcyc_post', 'trial_pow_post', '-v7.3')
%%

load Z:\HCP\Enrique\Results\allsubs_pow-dur_epis_LeftHand.mat
trial_pow_pre_Left = trial_pow_pre;
trial_pow_post_Left = trial_pow_post;
trial_durcyc_pre_Left = trial_durcyc_pre;
trial_durcyc_post_Left = trial_durcyc_post;
load Z:\HCP\Enrique\Results\allsubs_pow-dur_epis_RightHand.mat
trial_pow_pre_Right = trial_pow_pre;
trial_pow_post_Right = trial_pow_post;
trial_durcyc_pre_Right = trial_durcyc_pre;
trial_durcyc_post_Right = trial_durcyc_post;

clear trial_durcyc_post trial_durcyc_pre trial_pow_post trial_pow_pre

%% Analysis
frex = exp(0.6:0.1:3.7); % 1.8 Hz to 40 Hz
nFrex = length(frex);
nVox = 1925; 

% Differences
Nsub = size(trial_durcyc_pre_Left, 1);

group_diff_pow = zeros(Nsub, nVox, 'single');
group_diff_cycles = zeros(Nsub, nVox, 'single');

for s = 1:Nsub
    s
    % Number of cycles to seconds
    frex_3d = reshape(frex, 1, 1, []);
    trial_durcyc_pre_Right{s} = trial_durcyc_pre_Right{s}  ./ repmat(frex_3d, [size(trial_durcyc_pre_Right{s},1), nVox, 1]);
    trial_durcyc_post_Right{s} = trial_durcyc_post_Right{s}  ./ repmat(frex_3d, [size(trial_durcyc_post_Right{s},1), nVox, 1]);
    trial_durcyc_pre_Left{s} = trial_durcyc_pre_Left{s}  ./ repmat(frex_3d, [size(trial_durcyc_pre_Left{s},1), nVox, 1]);
    trial_durcyc_post_Left{s} = trial_durcyc_post_Left{s}  ./ repmat(frex_3d, [size(trial_durcyc_post_Left{s},1), nVox, 1]);

    % Mean of trials
    tmp_powR = squeeze(mean(trial_pow_pre_Right{s}, 1));
    tmp_powL = squeeze(mean(trial_pow_pre_Left{s}, 1));
    tmp_cycR = squeeze(mean(trial_durcyc_pre_Right{s}, 1));
    tmp_cycL = squeeze(mean(trial_durcyc_pre_Left{s}, 1));

    % Zscore Normalize
    tmp_powR(tmp_powR == 0) = NaN;
    tmp_powR_z = normalize(tmp_powR', 2, "zscore")';

    tmp_powL(tmp_powL == 0) = NaN;
    tmp_powL_z = normalize(tmp_powL', 2, "zscore")';

    tmp_cycR(tmp_cycR == 0) = NaN;
    tmp_cycR_z = normalize(tmp_cycR', 2, "zscore")';

    tmp_cycL(tmp_cycL == 0) = NaN;
    tmp_cycL_z = normalize(tmp_cycL', 2, "zscore")';

    % Recover 0s instead of NaNs
    tmp_powR_z(isnan(tmp_powR_z)) = 0; tmp_powL_z(isnan(tmp_powL_z)) = 0;
    tmp_cycR_z(isnan(tmp_cycR_z)) = 0; tmp_cycL_z(isnan(tmp_cycL_z)) = 0;

    % Raw data (before normalization for plotting)
    tmp_powR(isnan(tmp_powR)) = 0; tmp_powL(isnan(tmp_powL)) = 0;
    tmp_cycR(isnan(tmp_cycR)) = 0; tmp_cycL(isnan(tmp_cycL)) = 0;

    sub_avg_pre_powR(s,:,:) = tmp_powR_z;
    sub_avg_pre_powL(s,:,:) = tmp_powL_z;
    sub_avg_pre_cyclesR(s,:,:) = tmp_cycR_z;
    sub_avg_pre_cyclesL(s,:,:) = tmp_cycL_z;

    sub_avg_pre_powR_raw(s,:,:) = tmp_powR;
    sub_avg_pre_powL_raw(s,:,:) = tmp_powL;
    sub_avg_pre_cyclesR_raw(s,:,:) = tmp_cycR;
    sub_avg_pre_cyclesL_raw(s,:,:) = tmp_cycL;
end

group_diff_pow    = sub_avg_pre_powL    - sub_avg_pre_powR;
group_diff_cycles = sub_avg_pre_cyclesL - sub_avg_pre_cyclesR;

mean_pre_powL = squeeze(mean(sub_avg_pre_powL_raw));
mean_pre_powR = squeeze(mean(sub_avg_pre_powR_raw));
mean_pre_cycL = squeeze(mean(sub_avg_pre_cyclesL_raw));
mean_pre_cycR = squeeze(mean(sub_avg_pre_cyclesR_raw));

%% Ttest with empiric data
% Power
mean_pow = squeeze(mean(group_diff_pow, 1));
sd_pow   = squeeze(std(group_diff_pow, 0, 1)) / sqrt(Nsub);
t_pow    = mean_pow ./ sd_pow;
t_pow(isnan(t_pow))       = 0;

% Cycles
mean_cycles = squeeze(mean(group_diff_cycles, 1));
sd_cycles   = squeeze(std(group_diff_cycles, 0, 1)) / sqrt(Nsub);
t_cycles    = mean_cycles ./ sd_cycles;
t_cycles(isnan(t_cycles)) = 0;

%% Permutation analysis
nPerms = 1000;

max_t_null_pow    = zeros(nPerms,1, 'single');
max_t_null_cycles = zeros(nPerms,1, 'single');

for p = 1:nPerms
    % Flip signs
    sign_flips = randi([0, 1], Nsub, 1) * 2 - 1;   
    perm_diff_pow    = group_diff_pow .* sign_flips;
    perm_diff_cycles = group_diff_cycles .* sign_flips;
    
    perm_mean_pow = squeeze(mean(perm_diff_pow, 1));
    perm_sd_pow   = squeeze(std(perm_diff_pow, 0, 1)) / sqrt(Nsub);
    perm_t_pow    = perm_mean_pow ./ perm_sd_pow;
    
    perm_mean_cycles = squeeze(mean(perm_diff_cycles, 1));
    perm_sd_cycles   = squeeze(std(perm_diff_cycles, 0, 1)) / sqrt(Nsub);
    perm_t_cycles    = perm_mean_cycles ./ perm_sd_cycles;
    
    % Keep t-max of this perm 
    max_t_null_pow(p)    = max(abs(perm_t_pow),[], 'all');
    max_t_null_cycles(p) = max(abs(perm_t_cycles), [], 'all');
    
end

alpha = 0.05;
sorted_null_pow    = sort(max_t_null_pow,1);
sorted_null_cycles = sort(max_t_null_cycles,1);

idx_95 = round((1 - alpha) * nPerms);

t_crit_pow    = sorted_null_pow(idx_95);
t_crit_cycles = sorted_null_cycles(idx_95);

sig_mask_pow    = abs(t_pow) >= t_crit_pow;
sig_mask_cycles = abs(t_cycles) >= t_crit_cycles;

t_pow_result = t_pow;
t_pow_result(~sig_mask_pow) = 0;
t_cycles_result = t_cycles;
t_cycles_result(~sig_mask_cycles) = 0;

alphaband = 16:21;
betaband = 24:29;
fband = betaband

sBOSC_nii(mean(t_pow_result(:,fband),2), ['ttest_power_beta']);
sBOSC_nii(mean(t_cycles_result(:,fband),2), ['ttest_cycles_beta']);

cfg = [];
cfg.colmap = slanCM('vik');
cfg.colim = [-10 10]
% cfg.interpmethod = 'linear'
sBOSC_sourcefig(sum(t_pow_result(:,fband),2),[],cfg)

[valmax, voxmax] = max(t_pow_result(:,fband));
[~,idxmax] = max(valmax);
[valmin, voxmin] = min(t_pow_result(:,fband));
[~,idxmin] = min(valmin);

figure('Color', 'w', 'Units', 'centimeters', 'Position', [5, 5, 15, 10]);
p1 = plot(frex, mean_pre_powL(voxmin(idxmin),:), 'Color', [0 0.4470 0.7410], 'LineWidth', 2.5); % MATLAB default blue
hold on;
p2 = plot(frex, mean_pre_powR(voxmin(idxmin),:), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2.5); % MATLAB default red
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Arial'; % Standard font for publications
ax.LineWidth = 1.2;
ax.TickDir = 'out';    % Ticks point outward
ax.Box = 'off';        % Removes the top and right axis lines
xlim([min(frex) max(frex)]);
xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Raw Power', 'FontSize', 14, 'FontWeight', 'bold');
lgd = legend([p1, p2], {'Contra', 'Ipsi'}, 'Location', 'northeast');
lgd.Box = 'off';
lgd.FontSize = 10;


figure('Color', 'w', 'Units', 'centimeters', 'Position', [5, 5, 15, 10]);
p1 = plot(frex, mean_pre_powR(voxmax(idxmax),:), 'Color', [0 0.4470 0.7410], 'LineWidth', 2.5); % MATLAB default blue
hold on;
p2 = plot(frex, mean_pre_powL(voxmax(idxmax),:), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2.5); % MATLAB default red
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Arial'; % Standard font for publications
ax.LineWidth = 1.2;
ax.TickDir = 'out';    % Ticks point outward
ax.Box = 'off';        % Removes the top and right axis lines
xlim([min(frex) max(frex)]);
xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Raw Power', 'FontSize', 14, 'FontWeight', 'bold');
lgd = legend([p1, p2], {'Contra', 'Ipsi'}, 'Location', 'northeast');
lgd.Box = 'off';
lgd.FontSize = 10;

% cycles
[valmax, voxmax] = max(t_cycles_result(:,fband));
[~,idxmax] = max(valmax);
[valmin, voxmin] = min(t_cycles_result(:,fband));
[~,idxmin] = min(valmin);

figure('Color', 'w', 'Units', 'centimeters', 'Position', [5, 5, 15, 10]);
p1 = plot(frex, mean_pre_cycL(voxmin(idxmin),:), 'Color', [0 0.4470 0.7410], 'LineWidth', 2.5); % MATLAB default blue
hold on;
p2 = plot(frex, mean_pre_cycR(voxmin(idxmin),:), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2.5); % MATLAB default red
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Arial'; % Standard font for publications
ax.LineWidth = 1.2;
ax.TickDir = 'out';    % Ticks point outward
ax.Box = 'off';        % Removes the top and right axis lines
xlim([min(frex) max(frex)]);
xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Raw Power', 'FontSize', 14, 'FontWeight', 'bold');
lgd = legend([p1, p2], {'Contra', 'Ipsi'}, 'Location', 'northeast');
lgd.Box = 'off';
lgd.FontSize = 10;

figure('Color', 'w', 'Units', 'centimeters', 'Position', [5, 5, 15, 10]);
p1 = plot(frex, mean_pre_cycR(voxmax(idxmax),:), 'Color', [0 0.4470 0.7410], 'LineWidth', 2.5); % MATLAB default blue
hold on;
p2 = plot(frex, mean_pre_cycL(voxmax(idxmax),:), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2.5); % MATLAB default red
ax = gca;
ax.FontSize = 12;
ax.FontName = 'Arial'; % Standard font for publications
ax.LineWidth = 1.2;
ax.TickDir = 'out';    % Ticks point outward
ax.Box = 'off';        % Removes the top and right axis lines
xlim([min(frex) max(frex)]);
xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Raw Power', 'FontSize', 14, 'FontWeight', 'bold');
lgd = legend([p1, p2], {'Contra', 'Ipsi'}, 'Location', 'northeast');
lgd.Box = 'off';
lgd.FontSize = 10;

