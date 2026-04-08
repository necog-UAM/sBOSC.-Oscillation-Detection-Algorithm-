function [source_forward, source_inverse, datasource] = sBOSC_beamforming (dataclean, cfg)

% -------------------------------------------------------- %
% 3.1. MRI normalization (output: mri_norm)
% 3.2. Head model
% 3.3. Forward model (output: source_forward)
% 3.4. Source reconstruction (output: source_inverse)
% 3.5 Reconstruction of source-level activity (output: datasource)
% -------------------------------------------------------- %

mri = cfg.mri;
figflag = cfg.figure;

%% 3.1. MRI normalization
% mri normalized (mrin) and transformation matrix (normtrans)
cfg            = [];
cfg.nonlinear  = 'no';
cfg.spmversion = 'spm8';
mrin           = ft_volumenormalise(cfg, mri);          % do you want to change the anatomical labels for the axes [Y, n]? Y (r,a,s,i)

% determine the affine source->template coordinate transformation (from fieldtrip-20180405)
normtrans = mrin.params.VG.mat * inv(mrin.params.Affine) * inv(mrin.params.VF.mat) * mrin.initial;
    
%% 3.2. Head model
% semi-realistic singleshell head model based on the implementation from Guido Nolte
cfg             = [];
segment         = ft_volumesegment(cfg,mri);        % extract brain surface
segment.anatomy = mri.anatomy;

cfg        = [];
cfg.method = 'singleshell';
vol        = ft_prepare_headmodel(cfg, segment);    % construct semi-realistic singleshell head model

%% 3.3. Forward model
% output: source_forward (contains mri, vol, grad, and the normalized grid and leadfields)
grad  = dataclean.grad;

% Load normalized template grid (10mm)
load standard_sourcemodel3d10mm % From fieldtrip templates
grid = sourcemodel;

% Adapt the normalized grid to each individual's brain space
posmni    = grid.pos;
pos = ft_warp_apply(inv(normtrans), grid.pos*10, 'homogenous') / 10;
grid.pos  = pos;
grid.unit = 'cm';

% Convert grad, vol and grid to common units (mm)
grad = ft_convert_units(grad, vol.unit);
grid = ft_convert_units(grid, vol.unit);

% Select only voxels within cortical mask (e.g. cerebellum is excluded)
% and corrected (inside the cortical surface projected with ft_sourceplot)
% created with select_corticalvox_aal and add all 3423 voxels.

% select the 3423 voxels inside the cube (load from sBOSC)
load source_template_10mm_3423.mat
grid.inside = source.inverse.inside;

% Compute leadfields for each grid's voxel
cfg             = [];
cfg.grid        = grid;
cfg.grad        = grad;
cfg.vol         = vol;
cfg.channel     = {'MEG'};
cfg.normalize   = 'no';
cfg.reducerank  = 2;
grid2           = ft_prepare_leadfield(cfg);

if strcmp(figflag,'yes')
    figure
    plot3 (grad.chanpos(:,1), grad.chanpos(:,2), grad.chanpos(:,3), '.','MarkerEdgeColor',[0.8 0 0],'MarkerSize',25), hold on
    plot3 (vol.bnd.pos(:,1), vol.bnd.pos(:,2), vol.bnd.pos(:,3), '.','MarkerEdgeColor',[0 0 0.8]), hold on
    plot3 (grid2.pos(grid2.inside,1), grid2.pos(grid2.inside,2), grid2.pos(grid2.inside,3), '+k')
end

% Save grad, vol, grid and mri in source_forward structure to be used later
source_forward      = [];
source_forward.vol  = vol;
source_forward.mri  = mrin;
source_forward.grad = grad;
source_forward.grid = grid2;

%% 3.4. Computation of beamforming weights
% output: source_inverse (contains beamforming weights in source.avg.filter)
cfg            = [];
cfg.covariance = 'yes';
datacov        = ft_timelockanalysis(cfg, dataclean);      % covariance matrix

% Compute spatial filters (in source.avg.filter)
cfg                   = [];
cfg.method            = 'lcmv';
cfg.grad              = source_forward.grad;
cfg.headmodel         = source_forward.vol;
cfg.grid              = source_forward.grid;
cfg.lcmv.fixedori     = 'yes';
cfg.lcmv.normalize    = 'no';
cfg.lcmv.projectnoise = 'yes'; 
cfg.lcmv.keepfilter   = 'yes';          % important: save filters to use them later
cfg.lcmv.lambda       = '10%';          % the higher the smoother
cfg.lcmv.reducerank   = 2;
source                = ft_sourceanalysis(cfg, datacov);

load standard_sourcemodel3d10mm 
source.avg.ori = {};
source.avg.mom = {};
source.avg.noisecov = {};
source.pos     = sourcemodel.pos;            % standard grid positions
source.inside  = source_forward.grid.inside;
source_inverse = source;

%% 3.5 Reconstruction of source-level activity
% output: datasource
time         = dataclean.time{1};
voxel_inside = find(source.inside==1);
Nvox         = length(voxel_inside);
Ntrial       = length(dataclean.trial);
datasource   = zeros(Nvox,length(time),Ntrial);
[~, selch] = match_str(dataclean.label, source_inverse.avg.label);

for i = 1:Nvox
    sourcefilt = source.avg.filter{voxel_inside(i)}(selch);
    for tr = 1:Ntrial
        datasource(i,:,tr) = sourcefilt * dataclean.trial{tr};
    end
end

% Organize source-reconstructed data in a new fieldtrip structure
for tr = 1:Ntrial
    dataclean.trial{tr}  = single(squeeze(datasource(:,:,tr)));
    dataclean.time{tr}   = time;
end
nsamples = cellfun(@length, dataclean.time);
dataclean.sampleinfo = [cumsum([1; nsamples(1:end-1)']) , cumsum(nsamples')];
dataclean.label      = {};
for i = 1:Nvox
    dataclean.label{i} = ['V' num2str(i)];
end

datasource = dataclean;

end
