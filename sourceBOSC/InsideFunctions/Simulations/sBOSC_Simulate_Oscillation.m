function [oscillation, event_memory] = sBOSC_Simulate_Oscillation(cfg)
%% sBOSC_simulate_oscillation

% Get cfg input
time = cfg.sim.length;
fsample = cfg.sim.fsample;
nTp  = round(time*fsample);
unique_vox = unique([cfg.events.voxel],'stable'); 
nVox = length(unique_vox);

% Initialize
simsignal    = zeros(nVox, nTp);
occupancy    = false(nVox, nTp);
event_memory = cell(1, length(cfg.events));

%% Event loop
for ev = 1:length(cfg.events)
    
    event = cfg.events(ev);
    vidx  = find(unique_vox == event.voxel);
    freq    = event.freq;
    cycles  = event.cycles;

    if isfield(event, 'shape') && ~isempty(event.shape)
        event_shape = event.shape;
    else
        event_shape = 'sine';
    end
    if isfield(event, 'coexist_with') && ~isempty(event.coexist_with) && event.coexist_with ~= 0
        coexist_idx = event.coexist_with;
    else
        coexist_idx = 0;
    end    

    % Length of event
    npnts = round((cycles / freq) * fsample);
    if npnts > nTp, npnts = nTp; end

    start_idx = [];

    % Coexist
    if coexist_idx ~= 0
        parent_event = coexist_idx;
        start_idx = event_memory{parent_event}(1);

        % Trim check
        if start_idx + npnts - 1 > nTp
            npnts = nTp - start_idx + 1;
        end

    else
        % No coexistent event
        max_start = nTp - npnts + 1;
        if max_start < 1
            warning('Event %d is longer than the length of the signal. Omitting.', ev);
            continue;
        end

        ct = 0;
        valid = false;

        % If event is centered
        if isfield(event, 'centered') && strcmpi(event.centered, 'yes')
            temp_start = round((nTp - npnts) / 2) + 1;
            temp_end   = temp_start + npnts - 1;

            if sum(occupancy(vidx, temp_start:temp_end)) == 0
                start_idx = temp_start;
                valid = true;
            else
                warning("This event couldn't be placed because the positions were occupied by other signal. Use coexistence 'yes' to place both signals at the same positions.")
            end
            % Force out of loop
            ct = 1000;
        end

        % Search for occupancy
        while ~valid && ct < 1000
            temp_start = randi([1, max_start]);
            temp_end   = temp_start + npnts - 1;

            % IF any other episode in this time
            if sum(occupancy(vidx, temp_start:temp_end)) == 0
                start_idx = temp_start;
                valid = true;
            end
            ct = ct + 1;
        end

        if ~valid
            warning("Signal is full of events. Event %d in voxel %d couldn't be placed after 1000 attempts. ", ev, vidx);
            continue;
        end
    end

    % Save in memory for future events
    end_idx = start_idx + npnts - 1;
    event_memory{ev} = [start_idx, end_idx];

    %% Oscillation generation

    freq_vec = repmat(freq, 1, npnts);
    phase = 2 * pi * cumsum(freq_vec) / fsample;
    
    switch lower(event_shape)
        case 'mu'
            wave = sin(phase) + 0.5 * sin(2*phase + pi/2);
            wave = wave ./ max(abs(wave));
        case 'sine'
            wave = sin(phase);
    end
    
    % Smooth, amplitude scaling and normalization
    taper = tukeywin(npnts, 0.5)';
    final_wave = wave .* taper;
    final_wave = (final_wave ./ rms(final_wave));
    
    % Sum signals in case of coexistence
    simsignal(vidx, start_idx:end_idx) = simsignal(vidx, start_idx:end_idx) + final_wave;
    
    % Update occupancy matrix
    occupancy(vidx, start_idx:end_idx) = true;
end

oscillation = simsignal;
end