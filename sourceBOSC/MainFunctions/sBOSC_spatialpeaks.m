function  [spatialpks, localpks] = sBOSC_spatialpeaks(powspctm, thshld)

%% sBOSC_spatialpeaks
% Finds local maxima in the power spectrum, thresholds them using the
% previously computed aperiodic percentile and find spatial maxima in the 3D brain volume.
%
% Description:
% Keeps only peaks which: exceed the aperiodic threshold, are local 
% maxima in the time-frequency domain, and are spatial maxima in the brain.
%
% Input Arguments:
% - powspctm : [4D matrix] Power spectrum [nTrials x nVox x nFrex x nTp].
% - thshld   : [2D matrix] Aperiodic power threshold [nVox x nFrex].
%
% Output Arguments:
% - spatialpks : [4D logical] Boolean array of 3d volume peaks.
% - localpks   : [4D logical] Boolean array of local peaks that exceed the
%                aperiodic threshold.

%% Step 4. Find local peaks in each vox-freq point across time
load source_template_10mm_3423.mat
source = source.inverse;

[nTrials, nVox, nFrex, nTp] = size(powspctm);

% Spatial configuration
inside = find(source.inside);
positions = source.pos(inside,:);
distance = 1.5; % 19 neighbours

% Find connectivity matrix of voxels
connmat = pdist2(positions, positions) <= distance;
maxconn = max(sum(connmat));
valid_vox = find(sum(connmat, 2) == maxconn);

% Adjust positions
positions_allvox = source.pos - min(source.pos) + 1;
xx = positions_allvox(inside, 1)';
yy = positions_allvox(inside, 2)';
zz = positions_allvox(inside, 3)';
dim = source.dim;
ind_voxels = sub2ind(dim, xx, yy, zz);

% initialize
spatialpks = false(nTrials, nVox, nFrex, nTp);
localpks   = false(nTrials, nVox, nFrex, nTp);
thshld_3d = reshape(thshld, [nVox, nFrex, 1]);

%% Trial loop
for trl = 1:nTrials
    pow_trl = reshape(powspctm(trl, :, :, :), [nVox, nFrex, nTp]);
    mask_ap_ths = pow_trl >= thshld_3d;

    %% Find local peaks & Apply Aperiodic Threshold
    for vidx = 1:length(valid_vox)
        v = valid_vox(vidx);
        pow_v = reshape(pow_trl(v, :, :), [nFrex, nTp]);

        % Local peaks
        [~, ~, locmx] = findlocalmax(pow_v, 1, []);

        % Aperiodic threshold
        mask_ap_v = reshape(mask_ap_ths(v, :, :), [nFrex, nTp]);

        % Keep peaks that are local maxima and exceed the aperiodic threshold
        localpks(trl, v, :, :) = locmx & mask_ap_v;
    end

    %% Select only tf points with spatial maxima

    % Frequency loop
    for f0 = 1:nFrex
        anypeak = reshape(localpks(trl, :, f0, :), [nVox, nTp]);
        if ~any(anypeak(:)), continue; end

        t_filt = find(any(anypeak, 1));
        pow_f0 = reshape(pow_trl(:, f0, t_filt), [nVox, length(t_filt)]);

        % Interpolate to 3D grid
        dinterp = pow_f0' * source.dtempl;

        % Time loop
        for idx = 1:length(t_filt)
            t0 = t_filt(idx);
            vol3d = reshape(dinterp(idx, :), [dim(1), dim(2), dim(3)]);

            % Find spatial peaks
            [~, ~, locmx] = findlocalmax(vol3d, 3, []);
            % Extract peaks for inside voxels
            spatial_pks_vox = locmx(ind_voxels);
            % Intersection of local and spatial peaks
            spatialpks(trl, :, f0, t0) = (spatial_pks_vox(:) & anypeak(:, t0))';
        end

    end
end

end





