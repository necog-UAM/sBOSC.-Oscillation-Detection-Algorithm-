function [connected_episodes, conepisocc] = sBOSC_connect_episodes(cfg, epis)
%% sBOSC_connect_episodes
% Connects adjacent oscillatory episodes that overlap in time and frequency.
%
% Description:
% This function iteratively merges previously identified oscillatory episodes 
% that are contiguous or overlapping. Upon merging, frequency and power are recalculated 
% using a duration-weighted average.
%
% Input Arguments:
% - cfg  : [struct]     Configuration structure with fields:
%       cfg.frex    : [vector] Frequencies of interest in Hz.
%       cfg.fsample : [scalar] Sampling rate in Hz.
%       cfg.time    : [scalar] Total number of time points per trial.
% - epis : [cell array] {nTrials, nVox} array of structs containing
%                       descriptives (freq, dur_sec, dur_cyc, timeps, power) 
%                       for each original episode.
%
% Output Arguments:
% - connected_episodes : [cell array] {nTrials, nVox} array of structs containing
%                                     descriptives (freq, dur_sec, dur_cyc, timeps, power) 
%                                     for the newly merged episodes.
% - conepisocc         : [4D logical] Boolean mask of connected episode occurrences 
%                                     with dimensions [nTrials x nVox x nFrex x nTp].

fsample = cfg.fsample;
frex = cfg.frex;
nFrex = length(cfg.frex);
[nTrials, nVox] = size(epis);
nTp = cfg.time;

% Initialize
connected_episodes = cell(size(epis));
conepisocc = false(nTrials, nVox, nFrex, nTp);

security_dist = ceil(fsample / min(cfg.frex)) + 5;
offset = double(cfg.time) + security_dist;

for v = 1:nVox
    %%  Cat trials
    cat_epis_v = [];
    for trl = 1:nTrials
        if isempty(epis{trl, v}); continue; end
        curr_ep = epis{trl, v};
        for ep = 1:length(curr_ep)
            curr_ep(ep).timeps = double(curr_ep(ep).timeps) + ((trl - 1) * offset);
            curr_ep(ep).trl_id = trl;
        end
        if isempty(cat_epis_v)
            cat_epis_v = curr_ep;
        else
            cat_epis_v = [cat_epis_v, curr_ep];
        end
    end

    if isempty(cat_epis_v); continue; end

    listepis = [1:length(cat_epis_v)];           % original episodes, 1st iteration
    checkepis = [];
    catepis_connected = cat_epis_v(1);
    Nit = 0;             % number of iterations
    while ~isempty(listepis)
        if Nit > 0
            listepis = checkepis;
            checkepis = [];
            cat_epis_v = catepis_connected;
            catepis_connected(1) = cat_epis_v(1);
        end

        ep1 = 1;
        ctep = 1;
        while ep1 <= length(listepis)
            cat_epis_v(ep1).timeps = single(cat_epis_v(ep1).timeps);
            t0 = cat_epis_v(ep1).timeps(1);
            cycntps = 1/cat_epis_v(ep1).freq*fsample;       % number of time points corresponding to 1 cycle
            tf = round(cat_epis_v(ep1).timeps(end) + 1/2*cycntps);     % + 1/2 missing cycle
            tep1 = [t0:tf];
            if ep1+1 <= length(cat_epis_v) && cat_epis_v(ep1).freq ~= 100
                for ep2 = ep1+1:length(cat_epis_v)
                    cat_epis_v(ep2).timeps = single(cat_epis_v(ep2).timeps);
                    t0 = round(cat_epis_v(ep2).timeps(1) - 1/2*cycntps);   % - 1/2 missing cycle
                    cycntps = 1/cat_epis_v(ep2).freq*fsample;   % number of time points corresponding to 1 cycle
                    tf = cat_epis_v(ep2).timeps(end);
                    tep2 = [t0:tf];
                    if isempty(intersect(tep1,tep2))
                        ct = 1;
                        break;
                    elseif ~isempty(intersect(tep1,tep2))
                        f10 = (cat_epis_v(ep1).freq-0.1*cat_epis_v(ep1).freq);
                        f1f = (cat_epis_v(ep1).freq+0.1*cat_epis_v(ep1).freq);
                        f1 = [round(f10,1):0.1:round(f1f,1)];
                        f20 = (cat_epis_v(ep2).freq-0.1*cat_epis_v(ep2).freq);
                        f2f = (cat_epis_v(ep2).freq+0.1*cat_epis_v(ep2).freq);
                        f2 = [round(f20,1):0.1:round(f2f,1)];
                        if isempty(intersect(f1,f2))
                            ct = 1;
                        elseif ~isempty(intersect(f1,f2))
                            catepis_connected(ctep).freq = (cat_epis_v(ep1).freq*length(cat_epis_v(ep1).timeps) + cat_epis_v(ep2).freq*length(cat_epis_v(ep2).timeps)) ./ (length(cat_epis_v(ep1).timeps) + length(cat_epis_v(ep2).timeps));
                            catepis_connected(ctep).timeps = unique([tep1 tep2]);
                            catepis_connected(ctep).dur_sec = length(catepis_connected(ctep).timeps)./fsample;
                            catepis_connected(ctep).dur_cyc = catepis_connected(ctep).freq .* catepis_connected(ctep).dur_sec;
                            catepis_connected(ctep).power = (cat_epis_v(ep1).power*length(cat_epis_v(ep1).timeps) + cat_epis_v(ep2).power*length(cat_epis_v(ep2).timeps)) ./ (length(cat_epis_v(ep1).timeps) + length(cat_epis_v(ep2).timeps));
                            catepis_connected(ctep).trl_id = cat_epis_v(ep1).trl_id;
                            if ep1 == 1 && (length(listepis)==1)      % only 1st episode is a candidate
                                catepis_connected(ep2)= [] ;
                            end
                            ct = 0;
                            break;
                        end
                    end
                end
            else
                ct=1;
            end
            if ct==0
                checkepis = [checkepis ep1];      % new episodes might have new connections; run the algorithm again
                ctep = ctep + 1;
                ep1 = ep1+1;
                cat_epis_v(ep2).freq = 100;           % mark episodes that have already been connected
            elseif ct==1
                catepis_connected(ctep) = cat_epis_v(ep1);
                ctep = ctep + 1;
                ep1 = ep1+1;
            end
        end
        Nit = Nit+1;
    end

    catepis_connected = catepis_connected([catepis_connected.freq] ~= 100);

    %% Recover trials
    for ep = 1:length(catepis_connected)
        trl = catepis_connected(ep).trl_id;

        % Substract the offset previously added
        catepis_connected(ep).timeps = catepis_connected(ep).timeps - ((trl - 1) * offset);

        % Recover trial time
        valid_t = catepis_connected(ep).timeps >= 1 & catepis_connected(ep).timeps <= cfg.time;
        catepis_connected(ep).timeps = catepis_connected(ep).timeps(valid_t);
        
        if isempty(catepis_connected(ep).timeps)
            continue;
        end

        temptrl = rmfield(catepis_connected(ep), 'trl_id');

        if isempty(connected_episodes{trl, v})
            connected_episodes{trl, v} = temptrl;
        else
            connected_episodes{trl, v}(end+1) = temptrl;
        end

        % Build conepisocc matrix
        [~, f_idx] = min(abs(cfg.frex - temptrl.freq));
        conepisocc(trl, v, f_idx, temptrl.timeps) = true;
    end
end

end