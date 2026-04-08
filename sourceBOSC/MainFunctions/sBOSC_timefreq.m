function [powspctm, thshld, frex, fsample] = sBOSC_timefreq(data, aperiodic, cfg)

%% sBOSC_timefreq
% Performs a center of the head bias correction, computes the power spectrum of data and aperiodic, and calculates the aperiodic
% threshold.

% Description:
% This function corrects the center of the head bias using the aperiodic root-mean-square.
% It uses fieldtrip ft_freqanalysis to extract the power
% spectrum of the source_reconstructed data and also of the aperiodic data,
% using the same parameters. After that it calculates a threshold based on
% a percent of the aperiodic power for each voxel and frequency, that will
% be used later to filter out datasource power spectrum points that do not
% exceed this value.

% Input Arguments:
% - data     : [struct] Source-reconstructed data in FieldTrip format.
% - aperiodic: [struct] Source-reconstructed aperiodic data in Fieldtrip
%              format
% - cfg      : [struct] Configuration structure with fields:
%
%        cfg.frex        : [vector] Frequencies of interest.
%
%        cfg.frex_window : [scalar] Number of cycles to calculate the
%                          temporal window width. Default 5;
%        cfg.apthshld   : [scalar] Percentile of aperiodic power spectrum
%                          to extract the threshold
%
% Output Arguments:
% - powspctm : [4D matrix] The power spectrum of the original data. 
%              Dimensions are [nTrials x nVoxels x nFrex x nTime]. 
%
% - thshld   : [vector] The calculated aperiodic power threshold for each 
%              frequency using the configurated percentile 
%              Dimensions are [nVox x nFrex]. 
%
% - frex     : [vector] The exact frequencies evaluated and returned by 
%              FieldTrip's ft_freqanalysis. 
%              Dimensions are [1 x nFrex].
%
% - fsample  : [scalar] The new sampling frequency of the data after 
%              downsampling (128 Hz).


% Defaults
if ~isfield(cfg, 'frex_window'), frex_window = 5; else; frex_window = cfg.frex_window; end
if ~isfield(cfg, 'apthshld'), apthshld = 95; else; apthshld = cfg.apthshld; end

% Get cfg options
Frex = cfg.frex;
nFrex = length(Frex);

%% Step 1. Prepare and correct signal
cfg = [];
cfg.begsample = 1;
cfg.endsample = size(aperiodic.trial{1},2);
data = ft_redefinetrial(cfg,data);

nTrials = length(data.trial);
for trl = 1:nTrials
    rms_aperiodic = rms(aperiodic.trial{trl}, 2);
    data.trial{trl} = data.trial{trl} ./ rms_aperiodic;
    aperiodic.trial{trl} = aperiodic.trial{trl} ./ rms_aperiodic;
end

% II. Downsample both signals
cfg = [];
cfg.resamplefs = 128;
data = ft_resampledata(cfg, data);
aperiodic = ft_resampledata(cfg, aperiodic);
fsample = data.fsample;

%% Step 2.Frequency analysis and threshold

nVox    = size(data.trial{1}, 1);
nTime   = length(data.time{1});

start_samples = (0:nTrials-1)' * nTime + 1;
end_samples   = start_samples + nTime - 1;
data.sampleinfo = [start_samples, end_samples];
aperiodic.sampleinfo = [start_samples, end_samples];

powspctm = zeros(nTrials, nVox, nFrex, nTime, 'single');

frex = zeros(1, nFrex);
thshld = zeros(nFrex, 1, 'single');

% I. Time-Frequency decomposition
for fx = 1:nFrex
  
    % I. Original signal
    cfg            = [];
    cfg.method     = 'mtmconvol';
    cfg.taper      = 'hanning';
    cfg.output     = 'pow';
    cfg.foi        = Frex(fx);
    cfg.toi        = 'all';
    cfg.t_ftimwin  = frex_window./Frex(fx);
    cfg.pad        = 'nextpow2';
    cfg.keeptrials = 'yes';
    freq           = ft_freqanalysis(cfg, data);
    
    powspctm(:,:,fx,:) = single(freq.powspctrm);
    frex(fx) = freq.freq;

    clear freq

    % II. Aperiodic signal
    freqaperiodic  = ft_freqanalysis(cfg, aperiodic);
    powspctmaperiodic = single(freqaperiodic.powspctrm);

    % Compute the 95% percentile to obtain threshold
    powspctmaperiodic = permute(powspctmaperiodic, [2, 1, 4, 3]);  
    powspctmaperiodic = reshape(powspctmaperiodic, [nVox, nTrials * nTime]);
    
    thshld(fx) = prctile(prctile(powspctmaperiodic, 95, 1), 95, 2);
  
    clear freqaperiodic powspctmaperiodic
end
warning('on', 'all');

thshld = repmat(thshld', nVox, 1);

clear data

end
