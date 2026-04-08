function [episodes, episocc] = sBOSC_episodes(spatialpks, powspctm, cfg)
%% sBOSC_episodes
% Extracts and filters valid oscillatory episodes from spatial peaks.
%
% Description:
% This function takes the spatial peaks matrix and groups them into
% connected components in the time-frequency domain. It performs a spatial 
% smooth across neighboring voxels, filters out clusters that do not meet the minimum duration threshold 
% (in cycles) and calculates descriptives (frequency, duration, power) for 
% the surviving episodes.
%
% Input Arguments:
% - spatialpks : [4D logical] Spatial peaks [nTrials x nVox x nFrex x nTp].
% - powspctm   : [4D matrix]  Power spectrum [nTrials x nVox x nFrex x nTp].
% - cfg        : [struct]     Configuration structure with fields:
%       cfg.frex       : [vector] Frequencies of interest in Hz.
%       cfg.fsample    : [scalar] Sampling rate in Hz.
%       cfg.min_cycles : [scalar] Minimum duration of an episode in cycles.
%
% Output Arguments:
% - episodes   : [cell array] {nTrials, nVoxIN} array of structs containing
%                             descriptives (freq, dur_sec, dur_cyc, timeps, power) 
%                             for each episode.
% - episocc    : [4D logical] Clean boolean mask of valid episode occurrences 
%                             with dimensions [nTrials x nVoxIN x nFrex x nTp].

Frex = cfg.frex;
fsample = cfg.fsample;
min_cycles = cfg.min_cycles;

[nTrials, nVox, nFrex, nTp] = size(spatialpks);

load source_template_10mm_3423.mat
source = source.inverse;
inside = find(source.inside);

% Spatial configuration

positions = source.pos(inside,:);
distance = 1.5; % 19 neighbours maximum

% Find connectivity matrix of voxels
connmat = pdist2(positions, positions) <= distance;

% Load 1925 matrix
load source_template_10mm_1925.mat
inside1925 = find(source.inverse.inside);

voxin = find(ismember(inside,inside1925)); % The 1925 voxels inside the 3924
nVoxIN = length(voxin);

timep_freq = fsample.*(1./Frex') ;    % time points for 1 cycles of each foi

conn = conndef(2,'maximal');

% Initialize
episodes = cell(nTrials,nVoxIN);
episocc = false(nTrials, nVoxIN, nFrex, nTp);

%% Trial loop

for trl=1:nTrials

    spatialpks_trl = squeeze(spatialpks(trl, :, :, :));
    pow_trl = squeeze(powspctm(trl, :, :, :));

    for vin=1:nVoxIN   % only the 1925 voxels inside the cortex to save time
        vidx = voxin(vin);
        neighvoxs = find(connmat(vidx,:)==1);

        tfsmooth = squeeze(logical(mean(spatialpks_trl(neighvoxs, :, :),1)));

        if ~any(tfsmooth(:)), continue; end

        CC = bwconncomp(tfsmooth,conn);
        L = labelmatrix(CC);

        ct=1;
        for i=1:CC.NumObjects
            if length(CC.PixelIdxList{i})>timep_freq(end)*1       % 1--- only evaluate clusters with duration > 1 cycles of highest freq
                temp = L==i;
                
                % Weighted mean to get cluster frequency
                epis_freq = sum(temp,2);
                [~,idx_freq] = max(epis_freq);
                epis_freq = Frex(idx_freq);                
                epis_tps = logical(sum(temp,1));
                epis_dur = sum(epis_tps);
                epis_dur_sec = epis_dur.*1./fsample;
                epis_dur_cyc = epis_dur_sec./(1./epis_freq);

                    [f,t]=ind2sub(size(L),CC.PixelIdxList{i});
                    v2 = repmat(vidx,[length(f),1]);
                if epis_dur_cyc >= min_cycles 

                    % Store descriptives
                    episodes{trl,vin}(ct).freq = single(epis_freq);
                    episodes{trl,vin}(ct).dur_sec = single(epis_dur_sec);
                    episodes{trl,vin}(ct).dur_cyc = single(epis_dur_cyc);
                    episodes{trl,vin}(ct).timeps = int32(find(epis_tps==1));
                    idx_pow = sub2ind(size(pow_trl), v2, f, t);
                    pow = pow_trl(idx_pow);
                    episodes{trl,vin}(ct).power = single(mean(pow(:)));       
                    ct=ct+1;
                end
            end
        end
    end

    % Reconstruct episocc matrix
    for vin = 1:nVoxIN
        for ep=1:length(episodes{trl,vin})
            fm = episodes{trl,vin}(ep).freq;
            fbin = dsearchn(Frex',fm);
            tpts = episodes{trl,vin}(ep).timeps;

            episocc(trl, vin, fbin, tpts) = true;
        end
    end

end
end



