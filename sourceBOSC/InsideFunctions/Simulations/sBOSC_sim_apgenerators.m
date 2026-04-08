% Aperiodic generators. 
% Randomly select voxels to be aperiodic generators that are at a safe
% distance from designed oscillatory sources.

function apgens = sBOSC_Sim_apgenerators(cfg)

rng(cfg.seed)

% Load source models
load source_template_10mm_3423.mat

inside = find(source.inverse.inside);
positions = source.inverse.pos(inside,:);
distance = 5;
Nvox = length(inside);

% Find connectivity matrix of voxels
connmat = pdist2(positions, positions) <= distance;

for it=1:cfg.repetitions
    gen = randperm(Nvox,cfg.numgens);
    while ~isempty(intersect(gen,find(sum(connmat(cfg.Nsources,:))))) % If near neighbor of any sine generator, repeat
            gen = randperm(Nvox,cfg.numgens);
    end
    apgens(it,:) = gen;
end

end