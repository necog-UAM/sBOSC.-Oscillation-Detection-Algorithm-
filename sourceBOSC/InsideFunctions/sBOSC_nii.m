function sBOSC_nii(datasource, filename)

%--------------------------------------------%
% Generates a .nii file and saves it in current folder.
%--------------------------------------------%

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
source_mni.avg.pow(source_mni.inside) = datasource;
source_mni.time = 1;

load standard_mri.mat

cfg = [];
cfg.parameter     = 'avg.pow';
cfg.downsample    = 2;
cfg.interpmethod  = 'nearest';
source_interp     = ft_sourceinterpolate(cfg, source_mni, mri);

cfg = [];
cfg.filetype  = 'nifti';
cfg.parameter = 'pow';
cfg.filename  = filename;
ft_sourcewrite(cfg, source_interp)

end