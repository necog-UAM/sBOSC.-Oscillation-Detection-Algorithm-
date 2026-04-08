function [pks,locs,locmx] = findlocalmax(d,ndim,thr)
%
% Find local maxima (only positive peaks) in 1D, 2D or 3D data
% Other options (e.g., findpeaks or islocalmax) are much slower
%
% Use as:
%   [pks,locs] = findlocalmax(d,ndim,thr)
%
% Input:
%   d: data
%   ndim: number of dimensions (1:line, 2:plane, 3:volume)
%   [for ndim == 1 the size of d could be either n-by-1 or n-by-obs, 
%   where obs is the number of independent observations
%   thr: threshold for local maxima
% 
% Output:
%   pks: local maxima
%   locs: indices of local maxima
%   locmx: logical matrix indexing local maxima
% 
% Almudena Capilla, UAM

if ~exist('thr','var')
    thr = [];
end

if ndim == 1
    if size(d,1) == 1 && size(d,2) > 1
        d = d';
    end
    dx  = diff(squeeze(d(1:end-1,:)),1,1);          
    dx1 = diff(squeeze(d(2:end,:)),1,1);
    sgn = dx.*dx1;
    sgn = [zeros(1,size(sgn,2)); sgn; zeros(1,size(sgn,2))];
    
    d2x = diff(squeeze(d),2,1);
    d2x = [zeros(1,size(d2x,2)); d2x; zeros(1,size(d2x,2))];
    
    if size(d,2) == 1
        pks   = d(sgn(:,1)<0 & d2x(:,1)<0, 1);        % only positive peaks (2nd derivative negative)
        locs  = find(sgn(:,1)<0 & d2x(:,1)<0);
        if ~isempty(thr)
            pks2 = pks(pks > thr);
            locs = locs(pks > thr);
            pks = pks2;
        end
        % locmx = sgn(:,1)<0 & d2x(:,1)<0;
        locmx = zeros(size(d));
        locmx(locs) = 1;
        
    else
        locs = {};
        pks  = {};
        locmx = zeros(size(d));
        for i = 1:size(d,2)
            locs{i} = find(sgn(:,i)<0 & d2x(:,i)<0);
            pks{i}  = d(sgn(:,i)<0 & d2x(:,i)<0, i);
            if ~isempty(thr)
                pks2{i} = pks{i}(pks{i} > thr);
                locs{i} = locs{i}(pks{i} > thr);
                pks = pks2;
            end
            if ~isempty(locs{i})
                locmx(locs{i},i) = 1;
            end
        end
        locmx = logical(locmx);
    end
    
elseif ndim == 2
    % local maxima in x axis
    df  = diff(squeeze(d(1:end-1,:)),1,1);          
    df1 = diff(squeeze(d(2:end,:)),1,1);
    sgn = df.*df1;
    sgn = [zeros(1,size(sgn,2)); sgn; zeros(1,size(sgn,2))];
    
    d2f = diff(squeeze(d),2,1);
    d2f = [zeros(1,size(d2f,2)); d2f; zeros(1,size(d2f,2))];
    
    xlocmx = zeros(size(d));
    for i = 1:size(d,2)
        locs = find(sgn(:,i)<0 & d2f(:,i)<0);
        if ~isempty(locs)
            xlocmx(locs,i) = 1;
        end
    end
    
    % local maxima in y axis
    df  = diff(squeeze(d(:,1:end-1)),1,2);          
    df1 = diff(squeeze(d(:,2:end)),1,2);
    sgn = df.*df1;
    sgn = [zeros(size(sgn,1),1) sgn zeros(size(sgn,1),1)];
    
    d2f = diff(squeeze(d),2,2);
    d2f = [zeros(size(d2f,1),1) d2f zeros(size(d2f,1),1)];
    
    ylocmx = zeros(size(d));
    for i = 1:size(d,2)
        locs = find(sgn(:,i)<0 & d2f(:,i)<0);
        if ~isempty(locs)
            ylocmx(locs,i) = 1;
        end
    end
    
    locmx = logical(xlocmx.*ylocmx);
    pks   = d(locmx);
    [i,j] = find(locmx==1);
    if ~isempty(thr)
        for p = 1:length(pks)
            if pks(p) < thr
                locmx(i,j) = 0;
            end
        end
        pks   = d(locmx);
        [i,j] = find(locmx==1);
    end
    locs  = [i,j];
    
elseif ndim == 3
    locmx = zeros(size(d));
    
    for z = 2:size(d,3)-1
        dz = d(:,:,z);
        
        % local maxima in x axis
        df  = diff(squeeze(dz(1:end-1,:)),1,1);
        df1 = diff(squeeze(dz(2:end,:)),1,1);
        sgn = df.*df1;
        sgn = [zeros(1,size(sgn,2)); sgn; zeros(1,size(sgn,2))];
        
        d2f = diff(squeeze(dz),2,1);
        d2f = [zeros(1,size(d2f,2)); d2f; zeros(1,size(d2f,2))];
        
        xlocmx = zeros(size(dz));
        for i = 1:size(d,2)
            locs = find(sgn(:,i)<0 & d2f(:,i)<0);
            if ~isempty(locs)
                xlocmx(locs,i) = 1;
            end
        end
        
        % local maxima in y axis
        df  = diff(squeeze(dz(:,1:end-1)),1,2);
        df1 = diff(squeeze(dz(:,2:end)),1,2);
        sgn = df.*df1;
        sgn = [zeros(size(sgn,1),1) sgn zeros(size(sgn,1),1)];
        
        d2f = diff(squeeze(dz),2,2);
        d2f = [zeros(size(d2f,1),1) d2f zeros(size(d2f,1),1)];
        
        ylocmx = zeros(size(dz));
        for i = 1:size(d,2)
            locs = find(sgn(:,i)<0 & d2f(:,i)<0);
            if ~isempty(locs)
                ylocmx(locs,i) = 1;
            end
        end
        
        zlocmx = logical(xlocmx.*ylocmx);
        pks   = dz(zlocmx);
        [i,j] = find(zlocmx==1);
        locs  = [i,j];
        
        for i = 1:length(pks)
            if d(locs(i,1),locs(i,2),z-1) < pks(i) && d(locs(i,1),locs(i,2),z+1) < pks(i)
                locmx(locs(i,1),locs(i,2),z) = 1;
            end
        end
    end
    
    locmx = logical(locmx);
    pks   = d(locmx);
    locs  = [];
    for z = 2:size(d,3)-1
        [i,j] = find(locmx(:,:,z)==1);
        locsz  = [i,j];
        
        if ~isempty(locsz)
            for i = 1:size(locsz,1)
                locs = [locs ; locsz(i,1),locsz(i,2),z];
            end
        end
    end
   
    locs2 = [];
    ct = 1;
    if ~isempty(thr)
        for p = 1:length(pks)
            if pks(p) >= thr
                locs2(ct,:) = locs(p,:);
                ct = ct+1;
            else
                locmx(locs(p,1),locs(p,2),locs(p,3)) = 0;
            end
        end
        pks  = d(locmx);
        locs = locs2;
    end
end


