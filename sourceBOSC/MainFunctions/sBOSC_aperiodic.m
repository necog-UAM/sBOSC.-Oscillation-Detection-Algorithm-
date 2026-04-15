function aperiodic_data = sBOSC_aperiodic(data, cfg)

%% sBOSC_aperiodic
% Estimates the aperiodic component using FOOOF and reconstructs it in the 
% time domain.

% Description:
% This function estimates the aperiodic (1/f) component of the signal by:
% 1- Segmenting the data into shorter chunks (continous) or using blocks of trials.
% 2- Calculating the power spectrum density (PSD) of that segments/trials.
% 3- Applying FOOOF algorithm to the PSD to obtain the aperiodic component.
% 4- Reconstructing the aperiodic component in the time domain via IFFT.

% Input Arguments:
% - data : [struct] Source-reconstructed data in FieldTrip format.
% - cfg  : [struct] Configuration structure with fields:
%
%        cfg.datatype       : [string] Type of data. Options: 'continuous' or 'trials'.
%
%        cfg.windowlength   : [scalar] (For 'continuous' data only). Length
%                             of the window to estimate aperiodic component in seconds.
%                             [Default: 'all'].
%
%        cfg.trialblocksize : [scalar] (For 'trials' only). Number of
%                             trials to average in the frequency domain to estimate the
%                             aperiodic component. [Default: 'all'].
%
%        cfg.overlap        : [scalar] (For 'continuous' only). Proportion of 
%                             overlap between consecutive windows (0 to 1). 
%                             [Default: 0.5].
%
%        cfg.fooof          : [struct] FOOOF user parameters. If not
%                             specified, Brainstorm defaults are used.
%                             - .freq_range       : [min max] Range to fit.
%                             - .aperiodic_mode   : 'fixed' (no knee) or 'knee'.
%                             - .max_peaks        : Max n. of peaks (Default: 3).
%                             - .peak_width_limits: [min max] Width limits (Hz).
%                             - .min_peak_height  : Min peak height (dB).
%                             - .peak_threshold   : Threshold for detection (std).
%
% Output:
% - aperiodic_data: time domain aperiodic component of the signal.

% -------------------------------------------------------- %
% 0. Input Handling
% -------------------------------------------------------- %
if nargin < 2; cfg = []; end

% FieldTrip Format
if ~isstruct(data) || ~isfield(data, 'trial') || ~isfield(data, 'fsample')
    error('sBOSC:InvalidInput', 'Input "data" must be a FieldTrip structure with .trial and .fsample fields.');
end

% -------------------------------------------------------- %
% 1. Defaults
% -------------------------------------------------------- %

% 1.1 Datatype.
if ~isfield(cfg, 'datatype')
    if length(data.trial) == 1
        cfg.datatype = 'continuous';
        warning('sBOSC:DefaultParam', 'cfg.datatype not specified. Defaulting to "continuous".');
    else
        cfg.datatype = 'trials';
        warning('sBOSC:DefaultParam', 'cfg.datatype not specified. Defaulting to "trials".');
    end
end

% 1.2 Window Length (If cfg.datatype = 'continuous')
if strcmp(cfg.datatype, 'continuous')
    if ~isfield(cfg, 'windowlength') || isempty(cfg.windowlength)
        cfg.windowlength = 'all'; 
    end
end

% 1.3 Trial Block Size (If cfg.datatype = 'trials')
if strcmp(cfg.datatype, 'trials')
    if ~isfield(cfg, 'trialblocksize') || isempty(cfg.trialblocksize)
        cfg.trialblocksize = 'all';
    end
end

% -------------------------------------------------------- %
% 2. Data preparation
% -------------------------------------------------------- %

fs = data.fsample;

% Initialize output
aperiodic_data = data;

% Build matrix of [Voxels x Time x Trials]
if strcmp(cfg.datatype, 'continuous')
    % Continuous
    datamatrix = data.trial{1};
elseif strcmp(cfg.datatype, 'trials')
    try
        datamatrix = cat(3, data.trial{:});
    catch
        error('sBOSC:DimError', 'Trials must have same length.');
    end
end

% Define
if strcmp(cfg.datatype, 'continuous')
    [nVox, nTime] = size(datamatrix);
    nTrials = 1;
else
    [nVox, nTime, nTrials] = size(datamatrix);
end

% -------------------------------------------------------- %
% 3. Switch continuous or trials
% -------------------------------------------------------- %

switch cfg.datatype

%% Case: continuous
    case 'continuous'

        if ischar(cfg.windowlength) && strcmp(cfg.windowlength, 'all')
            windowsamples = nTime;
            overlap       = 0;
            step_samples  = nTime;
        else
            windowsamples = round(cfg.windowlength * fs);
        if isfield(cfg, 'overlap')
            overlap = cfg.overlap;
        else
            overlap     = 0.5;
            cfg.overlap = overlap;
        end
            step_samples = round(windowsamples * (1 - overlap));
        end

        aperiodic_continuous = zeros(nVox, nTime);
        overlap_ct           = zeros(1, nTime);

    % Start index
    win_starts = 1 :step_samples:nTime;
    
    for w = 1:length(win_starts)
        idx_st  = win_starts(w);
        idx_end = idx_st + windowsamples - 1;
    
        % Absorb short tail into this window
        if w < length(win_starts)
            next_end = win_starts(w+1) + windowsamples - 1;
            if next_end > nTime
                idx_end = nTime;
            end
        end
    
        idx_end = min(idx_end, nTime);
    
        datablock       = datamatrix(:, idx_st:idx_end);
        aperiodic_block = aperiodic2time(datablock, fs, cfg);
    
        % hanning
        current_win_len = idx_end - idx_st + 1;
        crossfade_win = hanning(current_win_len)';

        aperiodic_block_han = aperiodic_block .* crossfade_win;
        aperiodic_continuous(:, idx_st:idx_end) = aperiodic_continuous(:, idx_st:idx_end) + aperiodic_block_han;
        overlap_ct(idx_st:idx_end) = overlap_ct(idx_st:idx_end) + crossfade_win;

        if idx_end == nTime; break; end
    end

        % Average overlapping regions
        aperiodic_continuous = aperiodic_continuous ./ overlap_ct;

        aperiodic_data.trial{1} = aperiodic_continuous;
        aperiodic_data.time{1}  = aperiodic_data.time{1}(1:nTime);

    %% Case: trials
    case 'trials'
        % Size of the blocks    
        if ischar(cfg.trialblocksize) && strcmp(cfg.trialblocksize, 'all')
            block_size = nTrials;
        else
            block_size = cfg.trialblocksize;
        end

        % Initialize aperiodic matrix
        aperiodic_trials = zeros(nVox, nTime, nTrials);        
        % Loop over blocks
        block_starts = 1 : block_size : nTrials;
        
        for b = 1:length(block_starts)
            tr_idx_st  = block_starts(b);
            tr_idx_end = min(tr_idx_st + block_size - 1, nTrials);
        
            % Absorb short tail into this block
            if b < length(block_starts)
                next_end = block_starts(b+1) + block_size - 1;
                if next_end > nTrials
                    tr_idx_end = nTrials;
                end
            end
        
            tr_idx          = tr_idx_st:tr_idx_end;
            datablock       = datamatrix(:, :, tr_idx);
            aperiodic_block = aperiodic2time(datablock, fs, cfg);
        
            aperiodic_trials(:, :, tr_idx) = aperiodic_block;
        
            if tr_idx_end == nTrials; break; end
        end
        
        % Trial structure 
        for t = 1:nTrials
            aperiodic_data.trial{t} = aperiodic_trials(:, :, t);
        end
end

% -------------------------------------------------------- %
% Main function: FFT->FOOOF->IFFT reconstruction
% -------------------------------------------------------- %

function aperiodic_timesignal = aperiodic2time(datainside, fs, cfg)

    [nVoxblock, nfft, nTrialsblock] = size(datainside);

    %% 1. FFT

    % Vectorize data, reshape to [nVox x nTrials x nTime]
    data_vector = permute(datainside, [1 3 2]); 
    data_vector = reshape(data_vector, [], nfft);

    data_fft = fft(data_vector, [] ,2);

    if mod(nfft, 2) == 0
        limit = nfft/2 + 1;
        is_even = true;
    else
        limit = (nfft+1)/2;
        is_even = false;
    end

    % Get amplitude
    amp_fft = abs(data_fft(:, 1:limit)) ./ nfft;
    frq = (0:limit-1) * (fs / nfft);
   
    % Single-Sided
    if is_even
        % Exclude DC and Nyquist
        amp_fft(:, 2:end-1) = 2 * amp_fft(:, 2:end-1);
    else
        % Exclude DC
        amp_fft(:, 2:end)   = 2 * amp_fft(:, 2:end);
    end

    % Calculate Power
    pow_fft = (amp_fft.^2) / 2;
    pow_fft(:, 1) = amp_fft(:, 1).^2; % DC
    if is_even; pow_fft(:, end) = amp_fft(:, end).^2; end % Nyquist

    % Return to original dimensions
    pow_fft_matrix = reshape(pow_fft, [nVoxblock, nTrialsblock, length(frq)]);
    data_fft_matrix = reshape(data_fft, nVoxblock, nTrialsblock, nfft);

    %% 2. FOOOF
   
    is_continuous = isfield(cfg, 'datatype') && strcmp(cfg.datatype, 'continuous');

    % Apply Welch to smooth the spectrum
    if is_continuous 
        data_wch = datainside(:, :, 1)'; 
        
        welch_win = round(min(2, (nfft/fs) / 2) * fs);
        if welch_win > nfft; welch_win = nfft; end
        
        % Apply welch
        [pxx, ~] = pwelch(data_wch, hanning(welch_win), [], nfft, fs);
        
        % Scale PSD back to power
        mean_pow_fft = (pxx') * (fs/nfft); 
        
    else 
        avg_pow = mean(pow_fft_matrix, 2); % Trial average
        mean_pow_fft = reshape(avg_pow, [nVoxblock, limit]); 
    end

    % Configure input for brainstorm wrapper
    TF = zeros(nVoxblock , 1, limit-1); % Do not work with DC 0. Input all freqs but first
    TF(:, 1, :) = double(mean_pow_fft(:, 2:end));

    % Get FOOOF wrapper configuration
    opts = get_fooof_options(frq, cfg);

    % FOOOF
    [~, fg]= process_fooof("FOOOF_matlab", TF, frq(2:end), opts, opts.hasOptimTools);

    % Extract peak fit
    peakfit = [fg.peak_fit]; 
    peakfit = reshape(peakfit, [length(fg(1).peak_fit), length(fg)])';

    low_freq_mask = frq(2:end) < 2; % To avoid fooof from finding infraslow <2 Hz peaks that understimate aperiodic
    peakfit(:, low_freq_mask) = 1;
    aperiodic_ratio = 1 ./ peakfit;

    %% 3. IFFT    

    aperiodic_timesignal = zeros(nVoxblock, nfft, nTrialsblock);

    for tr = 1:nTrialsblock
        raw_power = reshape(pow_fft_matrix(:, tr, 2:end), [nVoxblock, limit-1]);
        raw_aperiodic = raw_power .* aperiodic_ratio;

        % Take the DC 
        dc_component = reshape(pow_fft_matrix(:, tr, 1), [nVoxblock, 1]);
        recap_pow = [dc_component, raw_aperiodic]; 

        % Recover amplitude
        recap_ampl = recap_pow;
        if is_even
            recap_ampl(:,2:end-1) = recap_pow(:,2:end-1) ./ 2;
        else
            recap_ampl(:,2:end) = recap_pow(:,2:end) ./ 2;
        end

        recap_ampl = sqrt(recap_ampl) * nfft;

        dataphase_complex = data_fft_matrix(:, tr, 1:limit);
        dataphase_matrix = reshape(dataphase_complex, [nVoxblock, limit]);

        % Take the original phase (individually in case of trials)
        dataphase = angle(dataphase_matrix);
        
        positivefx = recap_ampl .* exp(1i * dataphase);

        % Mirror positive frequencies for IFFT (exlc. nyq)
        if is_even
            fullspectra = [ positivefx, conj(fliplr(positivefx(:,2:end-1)))];
        else
            fullspectra = [ positivefx, conj(fliplr(positivefx(:,2:end)))];
        end

        % IFFT
        aperiodic_timesignal(:,:,tr) = real(ifft(fullspectra, nfft, 2));
    end
end


%% ------------------------------------------------------------------------
%  FOOOF Configuration Helper
% ------------------------------------------------------------------------
    function opt = get_fooof_options(frequency_range, cfg)

        % Brainstorm dependency
        if exist('ft_hastoolbox', 'file')
            try
                ft_hastoolbox('brainstorm', 1);
            catch
            end
        end

        % Fix to Nyquist (If not added +0.5, it gets masked)
        frequency_range(end) = frequency_range(end) + 0.5;

        % 2. Get Brainstorm Defaults
        try
            defaultopts = getfield(process_fooof('GetDescription'), 'options');
            bs_avail = true;
        catch
            bs_avail = false;
            defaultopts.peakwidth.Value{1}     = [0.5, 12];
            defaultopts.maxpeaks.Value{1}      = 3;
            defaultopts.minpeakheight.Value{1} = 0.0;
            defaultopts.apermode.Value         = 'fixed';
            defaultopts.peaktype.Value         = 'gaussian';
            defaultopts.proxthresh.Value{1}    = 2;
            defaultopts.guessweight.Value      = 'none';
            defaultopts.sorttype.Value         = 'box';
            defaultopts.sortparam.Value        = 'param';
            defaultopts.sortbands.Value        = [];
        end

        % Use cfg.fooof
        if isfield(cfg, 'fooof'), user_opts = cfg.fooof; else, user_opts = []; end

        opt = [];

        % Helper function to pick: User over Default
        if exist('ft_getopt', 'file')
            val = @(name, def) ft_getopt(user_opts, name, def);
        else
            val = @(name, def) get_opt_manual(user_opts, name, def);
        end

        opt.freq_range          = val('freq_range', frequency_range([1 end]));
        opt.peak_width_limits   = val('peak_width_limits', defaultopts.peakwidth.Value{1});
        opt.max_peaks           = val('max_peaks',         defaultopts.maxpeaks.Value{1});

        % dB to B conversion
        if bs_avail
            raw_min_height = defaultopts.minpeakheight.Value{1};
        else
            raw_min_height = 0;
        end
        opt.min_peak_height     = val('min_peak_height',   raw_min_height/10);

        opt.aperiodic_mode      = val('aperiodic_mode',    defaultopts.apermode.Value);
        opt.peak_threshold      = val('peak_threshold',    2); % 2 std dev: parameter for interface simplification
        opt.return_spectrum     = val('return_spectrum',   1); % SPM/FT: set to 1
        opt.border_threshold    = val('border_threshold',  1); % 1 std dev: proximity to edge of spectrum, static in Python

        % Matlab-only options
        opt.power_line          = val('power_line',        'inf'); % for some reason it should be a string, if you don't want a notch, use 'inf'. Brainstorm's default is '60'
        opt.peak_type           = val('peak_type',         defaultopts.peaktype.Value);
        opt.proximity_threshold = val('proximity_threshold', defaultopts.proxthresh.Value{1});
        opt.guess_weight        = val('guess_weight',      defaultopts.guessweight.Value);
        opt.thresh_after        = val('thresh_after',      true); % Threshold after fitting always selected for Matlab (mirrors the Python FOOOF closest by removing peaks that do not satisfy a user's predetermined conditions)

        % Output options
        opt.sort_type  = defaultopts.sorttype.Value;
        opt.sort_param = defaultopts.sortparam.Value;
        opt.sort_bands = defaultopts.sortbands.Value;

        % Check Bounds
        if (any(opt.freq_range < 0) || opt.freq_range(1) >= opt.freq_range(2))
            error('FOOOF:Config', 'Invalid Frequency range in FOOOF options.');
        end

        % Optimization Tools Check
        opt.hasOptimTools = 0;
        if exist('fmincon', 'file')
            opt.hasOptimTools = 1;
            disp('Using constrained optimization, Guess Weight ignored.')
        end
    end

    function val = get_opt_manual(st, name, def)
        if isfield(st, name)
            val = st.(name);
        else
            val = def;
        end
    end


end


