%% Batch script to run sBOSC in OMEGA data. 

%% ----------- 0. Prepare paths ---------------- %

% sBOSC path
p.sBOSC = 'Z:\Toolbox\sourceBOSC';
addpath(genpath(p.sBOSC))

% Fieldtrip path
addpath('Z:\Toolbox\fieldtrip-20230118')  % Fieldtrip path
ft_defaults

% OMEGA data path
p.data = 'Z:\OMEGA\OMEGA_data';
p.raw = 'Z:\OMEGA\OMEGA_raw';
p.scripts = 'Z:\OMEGA\Enrique\OMEGA-Repo\Scripts\OMEGA';
addpath(genpath(p.scripts))  % 

% Participant and session list
sub = [1 2 3 4 5 6 7 8 9 11 12 14 15 16 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 37 39 40 41 42 44 45 46 47 48 49 50 51 52 55 56 57 58 59 60 61 62 63 64 65 67 68 69 70 71 72 73 74 75 76 77 78 79 80 84 85 87 88 89 90 91 92 94 95 96 97 98 99 101 102 103 104 105 106 134 145 146 148 149 150 151 152 154 155 156 157 158 159 160 161 165 166 167 168 169 170 171 175 176 177 179 181 184 185 195 197 200 207 208 210 212]';
ses = [1 1 1 1 1 1 1 2 1  2  1  2  1  1  1  2  3  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  3  1  2  1  2  1  2  1  1  2  1  1  1  1  1  1  1  1  2  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1   1 ]';

subs = {} ; sess = {};
for s = 1:length(sub)
    if sub(s) < 10
        subs{s} = ['000' num2str(sub(s))];
    elseif sub(s) >= 10 & sub(s) < 100
        subs{s} = ['00' num2str(sub(s))];
    else
        subs{s} = ['0' num2str(sub(s))];
    end
    sess{s} = ['000' num2str(ses(s))];
end

Nsub  = length(sub);


%% ----------- 1. Reading, preprocessing and artifact correction ---------------- %

for s = 1:Nsub
    OMEGA_Preprocessing(subs{s}, sess{s}, p.raw, p.data);   
end

%% sourceBOSC
%% ----------------------- 2. Source reconstruction  ---------------------------- %

% MEG-MRI coregistration: mark fiducials (lpa, rpa and nasion; then quit) 
% and check result (if coregistration is not accurate, repeat procedure)
for s = 1:Nsub
    mri = ft_read_mri([p.raw '\sub-' subs{s} '\ses-' sess{s} '\anat\defaced_t1.nii']);
    cdfilepath([p.raw '\sub-' subs{s} '\ses-' sess{s} '\meg\'],'resting');      % 
    hsfile    = findfile('.pos');
    headshape = ft_read_headshape(hsfile);
    [mri,scp] = hscoreg([], mri, headshape);    % mark fiducials: lpa (l), rpa (r) and nasion (n), then quit (q)

    save([p.data '\sub-' subs{s} '\ses-' sess{s} '\mri_coreg.mat'], 'mri', 'scp');
end

% Compute forward model and beamforming weights
for s = 1:Nsub
    sprintf(['\n Processing Sub:' subs{s} '..\n'])
    load([p.data '\sub-' subs{s} '\ses-' sess{s} '\mri_coreg.mat'])
    load([p.data '\sub-' subs{s} '\ses-' sess{s} '\dataclean.mat'])   

    cfg = [];
    cfg.mri = mri;          % coregistered mri 
    cfg.figure = 'yes';     % to check that sensors, volyume, and grid are correctly coregistered
    [source_forward, source_inverse, datasource] = sBOSC_beamforming (dataclean, cfg);   % do you want to change the anatomical labels for the axes [Y, n]? Y (r,a,s,i)
     
    % Remove bad segments
    if exist([p.data '\sub-' subs{s} '\ses-' sess{s} '\badsegments.mat']) == 2     
        load([p.data '\sub-' subs{s} '\ses-' sess{s} '\badsegments.mat'])
        for b = 1:length(badsegments)
            [~, t1] = min(abs(datasource.time{1} - badsegments{b}(1)));
            [~, t2] = min(abs(datasource.time{1} - badsegments{b}(2)));
            datasource.time{1}(t1:t2) = [];
            datasource.trial{1}(:,t1:t2) = [];
            datasource.sampleinfo = [1 length(datasource.time{1})];
        end
    end

    % Use only 5-minute recordings for all participants
    datalength = 290;             
    if ~isempty(datalength)     
        if datasource.time{1}(end) > datalength      % >5 min
            [~, t1] = min(abs(datasource.time{1} - datalength));
            [~, t2] = min(abs(datasource.time{1} - datasource.time{1}(end)));
            datasource.time{1}(t1:t2) = [];
            datasource.trial{1}(:,t1:t2) = [];
            datasource.sampleinfo = [1 length(datasource.time{1})];
        end
    end

    % Resample data
    cfg = [];
    cfg.resamplefs = 256;
    datasource = ft_resampledata(cfg, datasource);

    save([p.data '\sub-' subs{s} '\ses-' sess{s} '\source_forward_10mm_3423.mat'], 'source_forward');
    save([p.data '\sub-' subs{s} '\ses-' sess{s} '\source_inverse_10mm_3423.mat'], 'source_inverse');
    save([p.data '\sub-' subs{s} '\ses-' sess{s} '\datasource_3423.mat'], 'datasource');
end

%% ----------- 3. Compute aperiodic to get thresholds ---------------- %
for s = 1:Nsub 
    sprintf(['\n Processing Sub:' subs{s} '..\n'])
    load([p.data '\sub-' subs{s} '\ses-' sess{s} '\datasource_3423.mat'])
    
        cfg = [];
        cfg.datatype   = 'continuous';
        cfg.windowlength = 20;
    evalc('aperiodic = sBOSC_aperiodic(datasource, cfg)');
    
    % Save
    save([p.data '\sub-' subs{s} '\ses-' sess{s} '\aperiodic.mat'], aperiodic);
    
    
%% ----------- 4. Power-spectrum and power threshold ---------------- %
        cfg = [];
        cfg.frex = exp(0.6:0.1:3.7);
        cfg.apthshld = 95;
    evalc('[powspctm, thshld, frex, fsample] = sBOSC_timefreq(datasource, aperiodic, cfg)');
    
    % Save
    save([p.data '\sub-' subs{s} '\ses-' sess{s} '\thshld.mat'], thshld);   
    
%% ----------- 5. Local and spatial peaks detection ---------------- %
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
    evalc('[conepisodes, conepisocc] = sBOSC_connect_episodes(cfg, episodes)');

    % Save (reduced size)
    mkdir([p.data '\sub-' subs{s} '\ses-' sess{s} '\episodes'])
    cd([p.data '\sub-' subs{s} '\ses-' sess{s} '\episodes'])
    for vin = 1:size(conepisodes,2)
        epis = conepisodes{vin};
        filename = sprintf('epis_v%d.mat', vin);
        save(filename, 'epis');
    end
   
end
