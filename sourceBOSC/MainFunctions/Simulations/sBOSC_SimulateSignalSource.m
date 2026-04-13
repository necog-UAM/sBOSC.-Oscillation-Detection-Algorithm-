function simsignal = sBOSC_SimulateSignalSource(cfg)
%% sBOSC_SimulteSignal_Source
% Simulates MEG signal in brain sources and projects it to MEG sensors.

% Description:
% This function simulates:
% 1- Aperiodic (1/f) signal at selected voxels inside the brain.
% 2- Oscillatory events at selected voxels inside the brain.
% Both components are scaled, projected to MEG sensors, and white noise is added.
%
% Input Arguments:
%
% cfg : struct
%   Configuration structure with the following fields:
%
%   General Settings:
%   length       : [numeric] Signal duration in seconds. (Default: 60)
%   fsample      : [numeric] Sampling frequency in Hz. (Default: 512)
%   apgenerators : [string or array] Defines where the 1/f background is 
%                  generated. Can be 'ROI' (groups voxels by AAL atlas, one random voxel from the ROI generates the aperiodic) 
%                  or an array of specific voxel indices (1 to 3423). (Default: 'ROI')
%   figures      : [string] 'yes' or 'no' to display simulation plots. 
%                  (Default: 'no')
%
%   Events Configuration:
%      events       : [array of structs] Defines the oscillatory events. 
%                     Each element must contain:
%       - .voxel        : [numeric] Index of the source voxel for the event.
%       - .freq         : [numeric] Frequency of the oscillation in Hz.
%       - .cycles       : [numeric] Number of cycles in the event.
%       - .snr          : [numeric] Target Signal-to-Noise Ratio (RMS signal / RMS noise).
%       - .snr_domain   : [string] Domain for the SNR calculation: SNR can
%                         be adjusted at the 'source' or 'channel'. If
%                         adjusted at the channel, SNR at source will be
%                         higher than input, since to achieve a certain SNR
%                         at the channels, the event at the source must be higher. 
%                         (Default: 'source').
%       - .shape        : [string] Waveform shape,'sine', 'mu'. 
%                         (Default: 'sine')
%       - .coexist_with : [numeric] Index of another event to temporally 
%                         overlap with. Set to 0 for independent timing.
%       - .centered     : [string] Whether the oscillation should be placed at the
%                         middle of the simulated signal. 
%                         (Default: 'no').
%                       
%
% Output Arguments:
% -----------------
% simsignal : struct
%   A FieldTrip data structure containing the simulated MEG data 
%   

%% Defaults

% Default signal_len of the signal is 60 seconds
if ~isfield(cfg, 'length'); signal_len = 60; else, signal_len = cfg.length; end
% Default sample rate is 512
if ~isfield(cfg,'fsample'); fsample = 512; else, fsample = cfg.fsample; end
% Default aperiodic generators are based on AAL ROIs
if ~isfield(cfg, 'apgenerators'); apgenerators = 'ROI'; else, apgenerators = cfg.apgenerators; end
% Default is to not generate figures
if ~isfield(cfg, 'figures')  figures = 'no'; else, figures = cfg.figures; end

%% Load a source template (from sBOSC)
load source_template_10mm_3423.mat

voxIN = find(source.inverse.inside);
nVox = length(voxIN);


%% IF aperiodic generators are based on AAL atlas, group voxels
if ischar(apgenerators) && strcmpi(apgenerators, 'ROI')
    roi_voxels = get_aal_rois(source, voxIN);
    nROIs = length(roi_voxels);
end

%% Get cfg input
simvox     = [cfg.events(:).voxel];
simfreq    = [cfg.events(:).freq];
simcycles  = [cfg.events(:).cycles];
simsnr     = [cfg.events(:).snr];
simcoexist = [cfg.events(:).coexist_with];
simshape   = string({cfg.events.shape});
nevents    = length([cfg.events]);
centered   = cell(1,nevents);
for nev = 1:nevents
    if isfield(cfg.events(nev), 'centered') && strcmpi(cfg.events(nev).centered, 'yes')
        centered{nev} = cfg.events(nev).centered;
    else
        centered{nev} = 'no'; % Default is not centered
    end
end
   
% Default simulation of SNR is based on source
if isfield(cfg.events, 'snr_domain'); simdomain = {cfg.events(:).snr_domain}; else, simdomain = repmat({'source'}, 1, nevents); end

%% 0. Get signal parameters

time  = 1/fsample:1/fsample:signal_len;
nTp  = length(time);

% Initialize brain matrix
aperiodic_vol = zeros(nVox, nTp);

%% 1. Generate aperiodic on brain sources

% Option: ROI.
if ischar(apgenerators) && strcmpi(apgenerators, 'ROI')
    for r = 1:nROIs
        vxs = roi_voxels{r};

        if ~isempty(vxs)
            roi_aperiodic = f_alpha_gaussian(nTp, 1, 1.2)';
            aperiodic_vol(vxs, :) = repmat(roi_aperiodic, length(vxs), 1);
        end
    end
    apgenlist = [roi_voxels{:}];
    napgens = length(apgenlist);

% Option: specific voxels
else
    for apvx = 1:length(apgenerators)
        aperiodic = f_alpha_gaussian(nTp, 1, 1.2)';
        aperiodic_vol(apgenerators(apvx), :) = aperiodic;
    end
        apgenlist = apgenerators;
        napgens = length(apgenerators);
end

    %--% Figure: Aperiodic data in sources
    if strcmp(figures,'yes')
        pinkcol = [255, 105, 180]/255;
        graycol = [0.8, 0.8, 0.8];
        fig = figure('Position', [0, 0, 1200, 900], 'Color', 'w');
        movegui(fig, 'center');
        subplot(4,2,[3 5])
        scatter3(source.inverse.pos(voxIN, 1), ...
                 source.inverse.pos(voxIN, 2), ...
                 source.inverse.pos(voxIN, 3), ...
                 20, graycol, 'filled', 'MarkerFaceAlpha', 0.3);
        hold on,
        scatter3(source.inverse.pos(voxIN(apgenlist), 1), ...
                 source.inverse.pos(voxIN(apgenlist), 2), ...
                 source.inverse.pos(voxIN(apgenlist), 3), ...
                 40, pinkcol, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1);
        axis equal; axis off; view(3);
        title('Simulated aperiodic sources')
        
        subplot(422)
        randap = randi(napgens);
        plot(time, aperiodic_vol(apgenlist(randap),:), 'Color', pinkcol, 'LineWidth',1)
        title('Aperiodic at the source')
    end

% Leadfield orientations
lfvx = zeros(length(source.forward.grid.label),nVox);
for vx = 1:nVox
    lfvx(:,vx) = source.forward.grid.leadfield{voxIN(vx)}(:,randi(3));
end

% Projection
aperiodicdata = lfvx*aperiodic_vol;
     
% Add gaussian noise on sensors
gausnoise = randn(size(aperiodicdata,1), size(aperiodicdata,2)) .* prctile(prctile(abs(aperiodicdata),50),50);
aperiodicdata = aperiodicdata + gausnoise;
aperiodicdata_norm = aperiodicdata ./ mean(rms(aperiodicdata,2),'omitmissing');

    %--% Figure: Aperiodic + gaussian noise  in sensors
    if strcmp(figures,'yes')
        subplot(424)
        plot(time, aperiodicdata_norm(randi(size(aperiodicdata_norm,1)),:),'Color', pinkcol, 'LineWidth',1)
        title('Projected aperiodic at the sensor + gaussian noise')
    end

%% 2. Generate oscillations on brain sources

% Initialize
oscillation_source = zeros(nVox, nTp);
oscillation_channel   = zeros(size(lfvx, 1), nTp);

cfg_osc = [];
cfg_osc.sim.fsample = fsample;
cfg_osc.sim.length = signal_len;
for ev=1:nevents
    cfg_osc.events(ev).voxel        = simvox(ev);
    cfg_osc.events(ev).freq         = simfreq(ev);
    cfg_osc.events(ev).cycles       = simcycles(ev);
    cfg_osc.events(ev).shape        = simshape(ev);
    cfg_osc.events(ev).coexist_with = simcoexist(ev);
    if strcmpi(centered{ev}, 'yes')
        cfg_osc.events(ev).centered = centered{ev};
    else
        cfg_osc.events(ev).centered = 'no'; % Default
    end
end

[oscillation, osc_pnts] = sBOSC_Simulate_Oscillation(cfg_osc);
unique_vox = unique(simvox, 'stable');

  %--% Figure: oscillation at source
    if strcmp(figures,'yes')
        subplot(426)
        plot(time, oscillation,'b', 'LineWidth',1)
        title('Oscillation at the source')
    end


%% Event loop (oscillations)
for ev = 1:nevents
    v_idx     = simvox(ev);
    ev_snr    = simsnr(ev);
    ev_domain = simdomain{ev};
    v_row = find(unique_vox == v_idx);

    % Sample points of the event
    start_idx = osc_pnts{ev}(1);
    end_idx   = osc_pnts{ev}(2);
    tps       = start_idx:end_idx;

    % This event
    ev_osc = zeros(1, nTp);
    ev_osc(tps) = oscillation(v_row, tps);

    % =========================================================
    % SNR can be defined at sources or at sensors
    % =========================================================
    switch lower(ev_domain)

        case 'source'
            v_ap_rms = rms(aperiodic_vol(v_idx, :));
            if v_ap_rms == 0
                v_ap_rms = mean(rms(aperiodic_vol(aperiodic_vol(:,1)~=0, :), 2));
            end

            % Scale based on SNR at the source
            scaled_osc = ev_osc * (v_ap_rms * ev_snr);

            % Store to project later
            oscillation_source(v_idx, :) = oscillation_source(v_idx, :) + scaled_osc;

        case 'channel'
            ev_vol = zeros(nVox, nTp);
            ev_vol(v_idx, :) = ev_osc;

            % Project to source just this event
            proj_ev = lfvx * ev_vol;

            % Channel with max projeciton
            [~, chidx] = max(rms(proj_ev(:, tps), 2));

            % RMS in this channel
            ev_rms = rms(proj_ev(chidx, tps));

            % Scale SNR in this channel
            if ev_rms > 0
                proj_ev_scaled = ev_snr .* proj_ev ./ ev_rms;
            else
                proj_ev_scaled = proj_ev;
            end

            % Store in channel space
            oscillation_channel = oscillation_channel + proj_ev_scaled;
    end
end

%% 3. Project to sensors

sim_sensors = lfvx * oscillation_source;
sim_sensors_norm = sim_sensors ./ mean(rms(aperiodicdata, 2), 'omitmissing');

% Sum aperiodic + gaussian + event
simdata = aperiodicdata_norm + sim_sensors_norm + oscillation_channel;

  %--% Figure: oscillation at source
    if strcmp(figures,'yes')
        [~, ch] = max(rms(sim_sensors_norm, 2));
        if ch >1; else, ch=chidx;end
        subplot(428)
        plot(time, simdata(ch,:),'b', 'LineWidth',1)
        title('Project all data to sensors')
    end

%% 4. Store data in a struct

simsignal = [];
simsignal.cfg = cfg;
simsignal.fsample = fsample;
simsignal.trial{1} = simdata;
simsignal.time{1} = time;
simsignal.label = source.forward.grad.label;
simsignal.label = ft_channelselection('M*', simsignal); % select magnmtrs

% Simulation parameters
simsignal.sim.events = cfg.events;
simsignal.sim.osctimepoints = osc_pnts;
simsignal.sim.apgenerators = apgenerators;
simsignal.sim.fsample = fsample;
simsignal.sim.oscillation = oscillation;

end

%% ========================================================================
%  Get AAL ROIS to generate aperiodic in each ROI
%  ========================================================================

function roi_voxels = get_aal_rois(source, voxIN)

ft_path = fileparts(which('ft_defaults'));
atlas_file = fullfile(ft_path, 'template', 'atlas', 'aal', 'ROI_MNI_V4.nii');
atlas = ft_read_atlas(atlas_file);

% Interpolate atlas to source
ft_warning('off', 'all');
cfg_interp = [];
cfg_interp.interpmethod = 'nearest';
cfg_interp.parameter    = 'tissue';
atlas_source = ft_sourceinterpolate(cfg_interp, atlas, source.inverse);
ft_warning('on', 'all');

% Group voxels by ROI
nROIs_total = length(atlas_source.tissuelabel);
roi_voxels = cell(nROIs_total, 1);

tissue = atlas_source.tissue(:);
inside = source.inverse.inside(:);


% Group inside voxels (3423) based on AAL ROI
for r = 1:nROIs_total
    vxs_absolute = find(tissue == r & inside);

    if ~isempty(vxs_absolute)
        [~, vxs_relative] = ismember(vxs_absolute, voxIN);
        randidx = randi(length(vxs_relative));
        roi_voxels{r} = vxs_relative(randidx); % Select one voxel for each ROI to generate the aperiodic
    end
end
end