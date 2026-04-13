function pBOSC_sourcefig(datasource, filename, varargin)

%--------------------------------------------%
% Generates a brain source image and saves it in current folder. Fieldtrip
% is needed. Can input color lims as cfg.colim.
%--------------------------------------------%

if nargin > 2
    optionalArg = varargin{1};
else
    optionalArg = 'default';
end

% Correct column orientation
dims = size(datasource);
if dims(2)>dims(1)
    datasource = datasource';
end

if dims(1) == 1925
    load source_template_10mm_1925.mat
    inside = source.inverse.inside;
    clear source
    
    load standard_sourcemodel3d10mm.mat 
    source_mni = sourcemodel;
    
elseif dims(1) == 3294
    load source_template_10mm_3294.mat
    inside = source.inverse.inside;
    clear source
    
    load standard_sourcemodel3d10mm.mat 
    source_mni = sourcemodel;
else
    error('Data input must have length of 1925 or 3294 voxels')
end

source_mni.inside = inside;

source_mni.avg.pow = nan(length(source_mni.inside), 1);
source_mni.avg.pow(inside) = datasource;
source_mni.time = 1;

load standard_mri.mat

cfg=[];
cfg.parameter  = 'avg.pow';
cfg.downsample = 2;
if strcmp(optionalArg, 'default')
    cfg.interpmethod = 'nearest';
elseif isfield(optionalArg,'interp')
    cfg.interpmethod = 'linear';
end
source_interp = ft_sourceinterpolate (cfg, source_mni, mri);

figure('WindowState','maximized','Color',[1 1 1]);
% figure

cfg               = [];
cfg.method        = 'surface';
cfg.funparameter  = 'pow';
cfg.maskparameter = cfg.funparameter;
if strcmp(optionalArg, 'default')
    cfg.funcolorlim   = 'auto';
elseif isfield(optionalArg,'colim')
    cfg.funcolorlim   = optionalArg.colim;
end
if strcmp(optionalArg, 'default')
    cfg.funcolormap   = 'auto';
elseif isfield(optionalArg,'colmap')
    cfg.funcolormap   = optionalArg.colmap;
end
cfg.projmethod    = 'nearest';
cfg.opacity       = 0.8;
cfg.figure        = 'gca';
cfg.camlight      = 'no';
cfg.colorbar      = 'yes';
cfg.surffile     = 'surface_pial_left.mat';
cfg.surfinflated  = 'surface_inflated_left_caret_white.mat';
subplot(2,2,1), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('left')
subplot(2,2,3), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('left')

cfg.surffile     = 'surface_pial_right.mat';
cfg.surfinflated  = 'surface_inflated_right_caret_white.mat';
subplot(2,2,2), ft_sourceplot(cfg,source_interp), view([90 0]),  camlight('right')
subplot(2,2,4), ft_sourceplot(cfg,source_interp), view([-90 0]), camlight('right')

if isfield(optionalArg,'savefig')
    print('-dtiff','-r300',filename);
end
end