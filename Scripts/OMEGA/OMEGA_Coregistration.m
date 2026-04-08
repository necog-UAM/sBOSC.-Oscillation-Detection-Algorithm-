function OMEGA_Coregistration (sub, ses, rawpath, dpath)

% sub: subect code, e.g., '0001'
% ses: session, e.g., '0001'
% rawpath: folder with raw data
% dpath: folder to store processed data

% -------------------------------------------------------- %
% 2.1. Coregistration of MEG-MRI spaces (output: mri_coreg)
% -------------------------------------------------------- %

%% 2.1. Coregistration of MEG-MRI spaces
% mri.transform is the transformation matrix to go from mri space to sensor space

mri = ft_read_mri([rawpath '/sub-' sub '/ses-' ses '/anat/defaced_t1.nii']);

cd([rawpath '\sub-' sub '\ses-' ses '\meg\'])
data = findfile('resting');
cd(data)                                  

hsfile    = findfile('.pos');
headshape = ft_read_headshape(hsfile);    
[mri,scp] = omega_coreg([], mri, headshape);    % mark fiducials: lpa (l), rpa (r) and nasion (n), then quit (q) 


save([dpath '\sub-' sub '\ses-' ses '\mri_coreg.mat'], 'mri', 'scp')

end