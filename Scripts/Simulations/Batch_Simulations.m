%% Batch script to run simulations. 

%% ----------- 0. Prepare paths ---------------- %

% sBOSC path
addpath(genpath('Z:\Toolbox\sourceBOSC'))

% Fieldtrip path
addpath('Z:\Toolbox\fieldtrip-20230118')  % Fieldtrip path
ft_defaults
p.data = 'Z:\OMEGA\Enrique\Simresults';

%% 1. Simulation parameters

% We select the following voxels based on its distance to sensors:
% 2924 Distance: 64.4301 mm
% 2711 Distance: 75.3871 mm
% 2469 Distance: 84.2666 mm

% Simulated voxels
Nsources = [2924, 2711, 2469];
% Simulated frequencies
Nfreqs =  [5 10 20];
% Simulated cycles
Ncycles = [3 10 20];
% Simulated SNR
Nsnr = [0.5 1 1.5];
% Number of repetitions
Nreps = 5;

% Voxels with aperiodic component. 
% Semi-randomly select voxels at a safe distance from osicllatory voxels
% We generate 5 random voxels for each repetition.
cfg = [];
cfg.numgens = 5;
cfg.repetitions = Nreps;
cfg.Nsources = Nsources;
cfg.seed = 211025; % seed for replication
apgens = sBOSC_sim_apgenerators(cfg);

%% Simulation loop
for vx = 1:length(Nsources) % Sources
    for fx = 1:length(Nfreqs) % Frequencies
        for cx = 1:length(Ncycles) % Cycles
            for snr = 1:length(Nsnr) % SNR
                for rep = 1:Nreps % repetitions
                    tic
                    fprintf([' \n Voxel: ' num2str(vx) '\n Frequency: ' num2str(fx) '\n Cycles: ' num2str(cx)  '\n SNR: ' num2str(snr) '\n Repetition: ' num2str(rep) '\n'])
                    cd(p.data)
                    mkdir(['Repetition_' num2str(rep)])
                    cd(['Repetition_' num2str(rep)])
                    cfg = [];
                    cfg.apgenerators = apgens(rep,:);
                    cfg.fsample  = 512;
                    cfg.length = 9;
                    cfg.figures = 'no';
                    cfg.events = [];
                    % Events
                    cfg.events(1).voxel  = Nsources(vx);
                    cfg.events(1).freq   = Nfreqs(fx);
                    cfg.events(1).cycles = Ncycles(cx);
                    cfg.events(1).snr = Nsnr(snr);
                    cfg.events(1).coexist_with = 0;
                    cfg.events(1).shape = 'sine';
                    cfg.events(1).snr_domain = 'channel';
                    cfg.events(1).centered = 'yes';
                simsignal = sBOSC_SimulateSignalSource(cfg);

                    %% Step 2: Source-reconstruction
                evalc('simsignal_source = sBOSC_SimulateBeamformer(simsignal)');
                    %% Step 3: Aperiodic
                    cfg = [];
                    cfg.datatype   = 'continuous';
                    cfg.windowlength = 3;
                evalc('sim_aperiodic = sBOSC_aperiodic(simsignal_source, cfg)');
                    %% Step 4: Powspctm and threshold
                    cfg = [];
                    cfg.frex = exp(0.6:0.1:3.7);
                    cfg.apthshld = 95;
                evalc('[powspctm, thshld, frex, fsample] = sBOSC_timefreq(simsignal_source, sim_aperiodic, cfg)');
                    %% Step 5: Spatial peaks
                    cfg = [];
                evalc('[spatialpks, localpks] = sBOSC_spatialpeaks(powspctm, thshld, cfg)');
                    %% Step 6: Construct episodes
                    cfg = [];
                    cfg.frex = frex;
                    cfg.fsample = fsample;
                    cfg.min_cycles = 1;
                evalc('[episodes, episocc] = sBOSC_episodes(spatialpks, powspctm, cfg)');
                    %% Step 7: Connect episodes
                    cfg = [];
                    cfg.fsample = fsample;
                    cfg.frex = frex;
                    cfg.time = size(spatialpks,4);
                evalc('[conepis, conepisocc] = sBOSC_connect_episodes(cfg, episodes)');
                    %% Step 8: Evaluate results
                model = sBOSC_Stats(conepisocc, simsignal_source, powspctm, fsample, frex, 0);

                    % Save model performance
                    save(['model_sbosc' num2str(vx), num2str(fx), num2str(cx), num2str(snr)], 'model');
                    reptime = toc;
                    fprintf([ ' \n Time of repetition: ' num2str( reptime / 60 ) ' min \n '])
                end
            end
        end
    end
end

