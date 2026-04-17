% Omega Natural frequencies

% sBOSC path
p.sBOSC = 'Z:\Toolbox\sourceBOSC';
addpath(genpath(p.sBOSC))

% Fieldtrip path
addpath('Z:\Toolbox\fieldtrip-20230118')  % Fieldtrip path
ft_defaults

% OMEGA data path
p.data = 'Z:\OMEGA\OMEGA_data';
p.raw = 'Z:\OMEGA\OMEGA_raw';

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

% Add jet omega mod color
addpath('Z:\OMEGA\Enrique\OMEGA-Repo\Scripts\OMEGA\Files')

% Frequency parameters
frex = exp(0.6:0.1:3.7); % 1.8 Hz to 40 Hz
Nfrex = length(frex);

load source_template_10mm_1925.mat

load([p.data '\oscillatory_results.mat'])

Nvoxin = size(freq_osc_time_cycles,2);
freq_osc_time_cycles = permute(freq_osc_time_cycles, [3, 2, 1]); %frex x vox x subs

%% Natural Frequencies of each participant
freqvoxs = zeros(size(freq_osc_time_cycles,2), size(freq_osc_time_cycles,3));
for s = 1:size(freqvoxs,2)
    s
    single_frexvox = freq_osc_time_cycles(:,:,s);
    % Normalize across voxels
    single_frexvox(single_frexvox==0)=NaN;
    single_frexvox = normalize(single_frexvox,2, "zscore");
    for v = 1:Nvoxin
        if isnan(max(single_frexvox(:,v)))
            freqvoxs(v,s) = NaN;
        else
            [p] = findlocalmax(single_frexvox(:,v),1);
            if ~isempty(p)
            fc = find(single_frexvox(:,v) == max(p));
            else
             fc = find(single_frexvox(:,v) == max(single_frexvox(:,v))); % peak freq of the cluster with max Z-value
            end
            if length(fc) > 1 % if more than one value, assign mode over neighbouring voxels
                vnb = find(connmat(v,:));
                for nb=1:length(vnb) % how many neighbors
                    tmpfc = find(single_frexvox(:,vnb(nb))==max(single_frexvox(:,vnb(nb))));
                    fcneigh(nb) = tmpfc(1); %%
                end
                fc = mode(fcneigh);
            end
            % if freq_osc_time_cycles(fc,v,s)==0 % If 0 total oscillations
            %     freqvoxs(v,s) = NaN;
            % end
            freqvoxs(v,s) = frex(fc);
        end
    end
end
    
%% Mean of participants
nfmean = mean(freqvoxs,2, "omitmissing");
natfmeancil = prctile(freqvoxs,2.5,2);          % 95% CI interval
natfmeanciu = prctile(freqvoxs,97.5,2);
natfmeanci = [natfmeancil natfmeanciu];
natfreqstd = std(freqvoxs,[],2,"omitmissing");
sBOSC_nii(nfmean, ['nf_mean']);
sBOSC_nii(natfreqstd, 'nf_std');

%% Gaussian fit

f1 = 0;
f2 = 3.7;
d  = 0.2;

foi1 = exp(f1:d:f2);
foi2 = exp(f1:d/10:f2);

for i = 1:length(foi2)
    c0 = foi2(i);
    rectp(i,:) = rectangularPulse(c0-c0/4 ,c0+c0/4 ,foi2);      % to isolate individual peaks later
end

par = parpool(4);                  % parallel pool
par.IdleTimeout = inf;

Nvox = size(freqvoxs,1);
Nboot = 500;
natf.A       = NaN(Nvox,Nboot);
natf.mu      = NaN(Nvox,Nboot);
natf.sigma   = NaN(Nvox,Nboot);
natf.rsq     = NaN(Nvox,Nboot);
natf.randsub = NaN(Nboot,Nsub);

% bootstrap confidence interval; based on bootci -> ci = bootci(Nboot,{@bootfun,freqvoxs},'type','per')
rng('shuffle')
for b = 1:Nboot
    disp(['Bootstrapping ' num2str(b) '/' num2str(Nboot)])
    randsub = randi(Nsub,Nsub);
    randsub = randsub(1,:);  % sample with replacement Nsub
    for v = 1:Nvox
        fv = freqvoxs(v,randsub);
        [counts,~] = histcounts(fv(:),foi1);
        centers = exp(log(foi1) + d/2);
        % figure(1), clf, plot(centers(1:end-1), counts,'--')        
        
        counts  = interp1(centers(1:end-1),counts,foi2,'pchip');
        centers = foi2;

        [pks,locs] = findpeaks(counts);
        if isempty(pks)
            [pks, locs] = max(counts);
        end
        [~,ii] = sort(pks,'descend');
        if length(ii) >= 3
            c0s = centers(locs(ii(1:3))); % fit a maximum of 3 peaks to speed up
        else
            c0s = centers(locs);
        end

        % figure(1), hold on, plot(centers, counts)
        % gaussian fit of all candidate peak frequencies
        % [A,mu,sigma,rsq] = gausfitc0_serial(c0s,counts,centers,rectp,foi2);
        [A,mu,sigma,rsq] = gausfitc0(c0s,counts,centers,rectp,foi2);
        
        % i = find(A == max(A));% identify the one with the highest goodness of fit 
        i = find(rsq == max(rsq));% identify the one with the highest goodness of fit
        i=i(1);
        natf.A(v,b)       = A(i);
        natf.mu(v,b)      = mu(i);
        natf.sigma(v,b)   = sigma(i);
        natf.rsq(v,b)     = rsq(i);
        natf.randsub(b,:) = randsub;
        gausf = A(i) * exp(-(centers-mu(i)).^2 /(2*sigma(i).^2));
        % plot(centers,squeeze(gausf),'k', 'LineWidth',2)
        % set(gca,'XLim',[0 35],'YLim',[0 80])
        % hold off
        % pause
    end
end

natfmu      = median(natf.mu,2,"omitmissing");
natfstd = std(natf.mu,[],2, "omitmissing");
natfci(:,1) = prctile(natf.mu,2.5,2);          % 95% CI interval
natfci(:,2) = prctile(natf.mu,97.5,2);

figure, plot(natfci(:,1)), hold on, plot(natfci(:,2))


cd('Z:\OMEGA\Enrique\Results')
save([p.data '\natfreq_bootstrap'], 'natf', 'natfmu', 'natfci', 'natfstd')


natfmu_btstrp = natfmu; % rename
natf_btstrp = natf;

% Figures
sBOSC_nii(natfmu_btstrp, 'natf_orig_bootstrap')

addpath('Z:\Toolbox\fieldtrip-20230118') % Path to Fieldtrip
colcfg = [];
colcfg.colim   = [0.7 3.4];  
colcfg.colmap   = 'jet_omega_mod';
sBOSC_sourcefig(log(natfmu_btstrp), 'Natfreq',colcfg)

%% Compare natfreq 
load natfreq_Nk25.mat % orighinal
natfmu_Almu      = median(natf.mu,2,"omitmissing");
natf_Almu = natf;

% Source plot
colcfg = [];
colcfg.colim   = [0.7 3.4];  
colcfg.colmap   = 'jet_omega_mod';
sBOSC_sourcefig(log(natfmu_Almu), 'test',colcfg)

% Correlacion

[rho,  pval] = corr(natfmu_btstrp, natfmu_Almu) %0.5404. 
[rho, pval] = corr(natfmu_aal, natfmu_aalMU) % 0.5083
[rho, pval] = corr(roinf, roinfMU) % 0.65 

figure, plot(natfmu_btstrp)
hold on, plot(natfmu_Almu)
legend({'pBOSC', 'Almu'})


%% Natural frequencies ROIs
% Load roi str
load aal_voxel_label_10mm.mat
Nroi = length(aal_label_reduc);

% figure
for roi = 1:Nroi
    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    muvoxs = natf_btstrp.mu(voxs,:);
    
    Nk = 4;
    [idx,C,sumd,D] = kmeans(muvoxs,Nk,'Replicates',5,'MaxIter',200); 
    nvoxs = [sum(idx==1) sum(idx==2) sum(idx==3) sum(idx==4) ];         
  Nk = Nk - sum(nvoxs==1);                                                                                    
    [idx,C,sumd,D] = kmeans(muvoxs,Nk,'Replicates',5,'MaxIter',200);    
    nvoxs = [sum(idx==1) sum(idx==2) sum(idx==3) sum(idx==4) ];                                
    dominclust = find(nvoxs == max(nvoxs));                              
    centroidvox = find(D(:,dominclust)==min(D(:,dominclust)));          
    represvox(roi) = voxs(centroidvox(1));                               
     % hold on, histogram(C'),title(aal_label_reduc{roi}),pause               
end

natfmu_aal      = median(natf_btstrp.mu(represvox,:),2);
natfci_aal(:,1) = prctile(natf_btstrp.mu(represvox,:),2.5,2);       % 95% CI interval
natfci_aal(:,2) = prctile(natf_btstrp.mu(represvox,:),97.5,2);

% load natfreq_aal
% save natfreq_aal natfmu_aal natfci_aal represvox aal_label_reduc natf

figure, plot(natfci_aal)

% Roi aissign
roinf = zeros(1925,1);
for roi = 1:Nroi
    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    roinf(voxs) = natfmu_aal(roi);
end

% Rois with pbosc
sBOSC_nii(roinf, 'Roi_NatfreqKmeans_orig')

% Hist allvox per ROI with pBOSC
figure,
sgtitle('All vox')
for rr = 1:Nroi
    subplot(4,10,rr)
% boxplot(natf_aal(rr,:))
inds = find(label_inside_aal_reduc==rr);
voxs = voxel_inside_aal(inds);
histogram(natf_btstrp.mu(voxs,:),frex)
xlim([0 35]), ylims(rr,:) = ylim();
title(aal_label_reduc{rr},'FontSize',8),
end

% Hist representative vox per ROI pBOSC
figure,
sgtitle('Representative vox')
for rr = 1:Nroi
    subplot(4,10,rr)
% boxplot(natf_aal(rr,:))
histogram(natf_btstrp.mu(represvox(rr),:),frex)
xlim([0 35])
title(aal_label_reduc{rr},'FontSize',8),
end

%% Natural frequencies ROIs Almu
% Load roi str
Nroi = length(aal_label_reduc);

% figure
for roi = 1:Nroi
    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    muvoxs = natf_Almu.mu(voxs,:);
    
    Nk = 4;
    [idx,C,sumd,D] = kmeans(muvoxs,Nk,'Replicates',5,'MaxIter',200);   
    nvoxs = [sum(idx==1) sum(idx==2) sum(idx==3) sum(idx==4) ];          
  Nk = Nk - sum(nvoxs==1);                                                                                    
    [idx,C,sumd,D] = kmeans(muvoxs,Nk,'Replicates',5,'MaxIter',200);    
    nvoxs = [sum(idx==1) sum(idx==2) sum(idx==3) sum(idx==4) ];                          
    dominclust = find(nvoxs == max(nvoxs));                            
    centroidvox = find(D(:,dominclust)==min(D(:,dominclust)));          
    represvox(roi) = voxs(centroidvox(1));                             
     % hold on, histogram(C'),title(aal_label_reduc{roi}),pause            
end

natfmu_aalMU      = median(natf_Almu.mu(represvox,:),2);
natfci_aalMU(:,1) = prctile(natf_Almu.mu(represvox,:),2.5,2);       % 95% CI interval
natfci_aalMU(:,2) = prctile(natf_Almu.mu(represvox,:),97.5,2);

figure, plot(natfci_aal)

% Roi aissign
roinf = zeros(1925,1);
for roi = 1:Nroi
    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    roinf(voxs) = natfmu_aalMU(roi);
end

% Hist allvox per ROI with pBOSC
for rr = 1:Nroi
    subplot(4,10,rr)
    hold on
% boxplot(natf_aal(rr,:))
inds = find(label_inside_aal_reduc==rr);
voxs = voxel_inside_aal(inds);
histogram(natf_Almu.mu(voxs,:),frex)
xlim([0 35]), ylims(rr,:) = ylim();
title(aal_label_reduc{rr},'FontSize',8),
end

% Hist representative vox per ROI
figure,
sgtitle('Representative vox')
for rr = 1:Nroi
    subplot(4,10,rr)
% boxplot(natf_aal(rr,:))
histogram(natf_Almu.mu(represvox(rr),:),frex)
xlim([0 35])
title(aal_label_reduc{rr},'FontSize',8),
end

%% Most differences  by ROI

diffs = abs(natfmu_aal - natfmu_aalMU);

data = [log(natfmu_aal) log(natfmu_aalMU)];
new_data = [natfmu_aal natfmu_aalMU]

labels = arrayfun(@(x) sprintf('%.1f', x), new_data(:), 'UniformOutput', false);
labels = reshape(labels, size(data));

figure('Position',[0 0 350 1000])
imagesc(data)
cmap = colormap("jet_omega_mod");
colorbar

% Normalize data 
dataNorm = (data - min(data(:))) / (max(data(:)) - min(data(:)));
[nRows, nCols] = size(data);

for row = 1:nRows
    for col = 1:nCols
        val = dataNorm(row, col);
        cmapIndex = max(1, round(val * (size(cmap,1)-1)) + 1);
        bgColor = cmap(cmapIndex, :);

        brightness = 0.299*bgColor(1) + 0.587*bgColor(2) + 0.114*bgColor(3);

        if brightness < 0.5
            textColor = 'white';
        else
            textColor = 'black';
        end

        text(col, row, labels{row, col}, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 12, ...
            'Color', textColor)
    end
end
set(gca, 'YTick', 1:nRows, 'YTickLabel', aal_label_reduc)
set(gca, 'XTick', 1:nCols)
axis tight

ax = gca;
ax.YTick = 1:size(data,1);              
ax.YTickLabel = aal_label_reduc;          
ax.TickLength = [0 0]; 
ax.YAxis.TickDirection = 'out';    

figure('Position',[0 0 150 1000])
imagesc(diffs)
cmap = flipud(gray);
colormap(cmap)
yticks([]), xticks([])

% Roi aissign
diffnf = zeros(1925,1);
for roi = 1:Nroi
    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    diffnf(voxs) = diffs(roi);
end
% diffnf=max(diffnf)-diffnf 
sBOSC_nii(diffnf, 'Diffs_ROIs_orig')

colcfg = [];
colcfg.colim   = [0 15];  
colcfg.colmap   = colormap('sky');
sBOSC_sourcefig(diffnf, 'test',colcfg)


















%% Revisar lo que hay a partir de aquí













% Load roi str
load(['Z:\OMEGA\OMEGA-NaturalFrequencies-main\mat_files\aal_voxel_label_10mm.mat'])
Nroi = length(aal_label_reduc);

f1 = 0;
f2 = 3.7;
d  = 0.2;

foi1 = exp(f1:d:f2);
foi2 = exp(f1:d/10:f2);

for i = 1:length(foi2)
    c0 = foi2(i);
    rectp(i,:) = rectangularPulse(c0-c0/4 ,c0+c0/4 ,foi2);      % to isolate individual peaks later
end

par = parpool(4);                  % parallel pool
par.IdleTimeout = inf;

natf_aal = zeros(Nroi, size(natf.mu,2));
for it = 1:size(natf.mu,2)
for roi = 1:Nroi
    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    fnatvox = natf.mu(voxs,it);
  
    [counts,~] = histcounts(fnatvox,foi1);
    centers = exp(log(foi1) + d/2);
    counts  = interp1(centers(1:end-1),counts,foi2,'pchip');
    centers = foi2;

    [pks,locs] = findpeaks(counts);
    [~,ii] = sort(pks,'descend');
    if length(ii) >= 3
        c0s = centers(locs(ii(1:3)));         % fit a maximum of 3 peaks to speed up
    else
        c0s = centers(locs);
    end

    % gaussian fit of all candidate peak frequencies
    [A,mu,sigma,rsq] = gausfitc0(c0s,counts,centers,rectp,foi2);
    % figure, hold on, plot(centers, counts)
    % plot(centers,squeeze(gausf),'k')

    i = find(rsq == max(rsq));      
    gausf = A(i) * exp(-(centers-mu(i)).^2 /(2*sigma(i).^2));
        % identify the one with the highest goodness of fit / amplitude
    natf_aal(roi,it)     = mu(i);
end
disp(['Iteration:' num2str(it) '/' num2str(size(natf.mu,2))])
end

roici(:,1) = prctile(natf_aal, 2.5,2)
roici(:,2) = prctile(natf_aal, 97.5,2)
roimu = median(natf_aal,2,"omitmissing");

roinf = zeros(1925,1);
for roi = 1:Nroi
    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    roinf(voxs) = roimu(roi);
end

save roinatfreqgauss3 roimu roici natf_aal roinf


Omega_nii(roinf, 'Roi_Natfreq')

% CIs
figure, plot(roici)
difci = roici(:,2) - roici(:,1);
roici = zeros(1925,1);
for roi = 1:Nroi
    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    roici(voxs) = difci(roi);
end
Omega_nii(roici, 'Roi_CI')



% Distribuciones
figure,
for rr = 1:Nroi
    subplot(4,10,rr)
% boxplot(natf_aal(rr,:))
histogram(natf_aal(rr,:))
title(aal_label_reduc{rr},'FontSize',8),
xlim([0 30])
end





















%% Natural frequencies (old)
freq_osc_time_cycles2 = freq_osc_time_secs;
freq_osc_time_cycles2(freq_osc_time_cycles2 == 0) = NaN;
all_natfreq_norm = normalize(freq_osc_time_cycles2,2,"zscore","robust");
avg_natfreq_norm = squeeze(mean(all_natfreq_norm,1,"omitmissing"));
nf= zeros(size(avg_natfreq_norm,1),1);

for v = 1:Nvoxin
    [pks,locs] = findpeaks(avg_natfreq_norm(v,:));
    if isempty(pks) | isinf(pks)
        nf(v) = NaN;
    else
        pks2 = pks(~isinf(pks));
        locs = locs(~isinf(pks));
        temp=frex(locs(find(pks2==max(pks2))));
        nf(v)=temp(1);
    end
end
% save natfeq.mat nf
Omega_nii(nf, 'Nf3cycmean')

% % Dominant frequencies
% avg_domfreq_norm = squeeze(mean(freq_osc_time_cycles,1));
% df = zeros(size(avg_domfreq_norm,1),1);
% 
% for v = 1:Nvoxin
%     [pks,locs] = findpeaks(avg_domfreq(v,:));
%     pks2 = pks(~isinf(pks));
%     locs = locs(~isinf(pks));
%     temp=frex(locs(find(pks2==max(pks2))));
%     df(v)=temp(1);
% end
% Omega_nii(nf, 'new_domfreqepis')

% Figure
colcfg.colim   = [0.7 3.4];  
colcfg.colmap   = 'jet_omega_mod';
Omega_sourcefig(log(natfmu), 'test',colcfg)


%% OLD
for roi = 1:Nroi
    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    muvoxs = natf.mu(voxs,:);
    muvoxs(isnan(muvoxs)) =0; %%%%%%%%%% quito nans porque kmeans no converge
    
    Nk = 4;
    [idx,C,sumd,D] = kmeans(muvoxs,Nk,'Replicates',5,'MaxIter',200);      % k-means 4 clusters max
    nvoxs = [sum(idx==1) sum(idx==2) sum(idx==3) sum(idx==4) ];           % number of voxels in each cluster                       
    Nk = Nk - sum(nvoxs==1);                                              % eliminate clusters with only 1 member and repeat clustering                                           
    [idx,C,sumd,D] = kmeans(muvoxs,Nk,'Replicates',5,'MaxIter',200);      % k-means adapting number of clusters 
    nvoxs = [sum(idx==1) sum(idx==2) sum(idx==3) sum(idx==4) ];           % number of voxels in each cluster                       

    dominclust = find(nvoxs == max(nvoxs));                               % dominant cluster
    centroidvox = find(D(:,dominclust)==min(D(:,dominclust)));            % voxel closer to centroid
    represvox(roi) = voxs(centroidvox(1));                                % representative voxel for this AAL region
     % histogram(C'),title(aal_label_reduc{roi}),pause                    % visualize pattern of natural frequencies for both clusters
end

natfmu_aal      = median(freqvoxs(represvox,:),2, "omitmissing");
natfci_aal(:,1) = prctile(freqvoxs(represvox,:),2.5,2);       % 95% CI interval
natfci_aal(:,2) = prctile(freqvoxs(represvox,:),97.5,2);
% [natfmu_aal natfci_aal]

roinf = zeros(1925,1);
for roi = 1:Nroi
    inds = find(label_inside_aal_reduc==roi);
    voxs = voxel_inside_aal(inds);
    roinf(voxs) = natfmu_aal(roi);
end

Omega_nii(roinf, 'Roi_Natfreq')


% Correlación
idx = find(roinf>0);

corr(natfmu(idx), roinf(idx))
figure,plot(natfmu(idx),'.')
hold on
plot(roinf(idx),'.')

