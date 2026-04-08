%% Step 2. Signal from sensors is reconstructed to source space using beamformers

function simsignal_source = sBOSC_SimulateBeamformer(simsignal)

% Load source models (from sBOSC)
load source_template_10mm_3423.mat
forward = source.forward;

grad  = forward.grad;
grid = forward.grid;
vol = forward.vol;

clear source

grad = ft_convert_units(grad, vol.unit);
grid = ft_convert_units(grid, vol.unit);

cfg             = [];
cfg.grid        = grid;
cfg.grad        = grad;
cfg.vol         = vol;
cfg.channel     = {'MEG'};
cfg.normalize   = 'no'; 
cfg.reducerank  = 2;
grid2           = ft_prepare_leadfield(cfg);

forward.grid = grid2;

cfg            = [];
cfg.covariance = 'yes';
datacov        = ft_timelockanalysis(cfg, simsignal); % covariance matrix

% Compute spatial filters (in inverse.avg.filter)
cfg                   = [];
cfg.method            = 'lcmv';
cfg.grad              = grad;
cfg.headmodel         = vol;
cfg.grid              = forward.grid;
cfg.lcmv.fixedori     = 'yes';
cfg.lcmv.normalize    = 'no'; 
cfg.lcmv.projectnoise = 'yes';
cfg.lcmv.keepfilter   = 'yes';          
cfg.lcmv.lambda       = '10%';        
cfg.lcmv.reducerank   = 2;
inverse                = ft_sourceanalysis(cfg, datacov);

load standard_sourcemodel3d10mm % from Fieldtrip
inverse.avg.ori = {};
inverse.avg.mom = {};
inverse.avg.noisecov = {};
inverse.pos     = sourcemodel.pos; % standard grid positions

%% 2.1 Reconstruction of source-level activity
time         = simsignal.time{1};
voxel_inside = find(inverse.inside);
Nvox         = length(voxel_inside);
datasource   = zeros(Nvox,length(time));
for i = 1:Nvox
    disp([num2str(i) ' / ' num2str(length(voxel_inside))])
    datasource(i,:) = inverse.avg.filter{voxel_inside(i)} * simsignal.trial{1};
end

% Organize source-reconstructed data in a new fieldtrip structure
simsignal.trial{1}   = single(datasource);
simsignal.time{1}    = time;
simsignal.label      = {};
simsignal.sampleinfo = [1 size(simsignal.trial{1},2)];

for i = 1:Nvox
    simsignal.label{i} = ['V' num2str(i)];
end
clear datasource

% Downsample signal
cfg = [];
cfg.resamplefs = 256;
[simsignal_source] = ft_resampledata(cfg, simsignal)

end
