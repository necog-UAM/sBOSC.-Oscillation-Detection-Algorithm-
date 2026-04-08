 %% Evaluate model

 function stats = sBOSC_Stats(predictedepisodes, simdata, powspctm, fsample, frex, figflag)
%
 if ~exist('figflag')
    figflag = 0;
 end

% Squeeze the trials dimension (only 1 trial in simulations)
predictedepisodes = squeeze(predictedepisodes);
powspctm = squeeze(powspctm);

simepisodes = false(size(powspctm));
dipolevox = simdata.sim.events.voxel;
dipolefrex = simdata.sim.events.freq;
dipolecyc = simdata.sim.events.cycles;
dwnsmpl = simdata.sim.fsample ./ fsample;
for simdips = 1:length(simdata.sim.events)
    ep_start = simdata.sim.osctimepoints{simdips}(1);
    ep_end   = simdata.sim.osctimepoints{simdips}(end);
    ep_dwn_start = max(1, round(ep_start / dwnsmpl));
    ep_down_end   = round(ep_end / dwnsmpl);
    dipolepnts = ep_dwn_start:ep_down_end;
    simepisodes(dipolevox(simdips), dsearchn(frex', dipolefrex(simdips)), dipolepnts) = 1;
end
load source_template_10mm_1925.mat
source1925.inside = find(source.inverse.inside);
load source_template_10mm_3423.mat
source3423.inside = find(source.inverse.inside);
voxin = find(ismember(source3423.inside,source1925.inside)); % The 1925 voxels inside the 3924

simulatedepisodes = simepisodes(voxin,:,:);
for simdips = 1:length(simdata.sim.events)
    dipolevoxIN(simdips) = find(voxin==dipolevox(simdips)); % Voxel with osillation in the 1925 inside voxs
end

%% Stats
diffepis = predictedepisodes + simulatedepisodes;
truepositive =  sum(sum(sum(diffepis==2))) ./ sum(sum(sum(simulatedepisodes))) * 100;
simfrx = find(sum(simulatedepisodes(dipolevoxIN,:,:),3));
simpnts = find(simulatedepisodes(dipolevoxIN,simfrx,:));

    % Figure: Model prediction
    if figflag
    cmp = colormap([1 1 1; [0 9 255]./255]);
    figure, imagesc(squeeze(predictedepisodes(dipolevoxIN,:,:))), axis xy, axis off
    colormap(cmp)
    figure, imagesc(squeeze(simulatedepisodes(dipolevoxIN,:,:))), axis xy, axis off
    colormap(cmp)
    end

% % Find the voxels with maximum number of the predicted episodes
[tmpval,leakvxidx] = sort(sum(predictedepisodes(:,simfrx,simpnts),3) / sum(simulatedepisodes(dipolevoxIN,simfrx,simpnts)),'descend');
leakvxidx = leakvxidx(tmpval>.9); % Voxels with at least 90% of episodes found

% Now keep the voxel with maximum power
    [~,maxvox] = max(sum(powspctm(leakvxidx,simfrx,simpnts),3));
    maxvox = leakvxidx(maxvox);

% Stats simulation vs prediction
diffepis = predictedepisodes - simulatedepisodes;
falsenegative = sum(sum(sum(diffepis==-1))) ./ sum(sum(sum(simulatedepisodes))) * 100;
falsepositive = sum(sum(sum(diffepis==1)))  ./ sum(sum(sum(simulatedepisodes==0)))  * 100;
truenegative =  sum(sum(sum(predictedepisodes==0))) ./  sum(sum(sum(simulatedepisodes==0))) * 100;
precision = truepositive / (truepositive + falsepositive) * 100;
recall = truepositive / (truepositive + falsenegative) * 100;
f1score = 2 * (precision*recall) ./ (precision + recall);

% Alternative voxel
altvxepis = squeeze(predictedepisodes(maxvox,simfrx,:)) + squeeze(simulatedepisodes(dipolevoxIN,simfrx,:));
truepositivealtvx =  sum(altvxepis==2) ./ sum( squeeze(simulatedepisodes(dipolevoxIN,simfrx,:))) * 100;
altvxepis = squeeze(predictedepisodes(maxvox,simfrx,:)) - squeeze(simulatedepisodes(dipolevoxIN,simfrx,:));
falsenegativealtvx =  sum(altvxepis==-1) ./ sum( squeeze(simulatedepisodes(dipolevoxIN,simfrx,:))) * 100;

truenegativealtvx =  sum(sum(sum(predictedepisodes==0))) ./  sum(sum(sum(simulatedepisodes==0))) * 100;
precisionaltvx = truepositivealtvx / (truepositivealtvx + falsepositive) * 100;
recallaltvx = truepositivealtvx / (truepositivealtvx + falsenegative) * 100;
f1scorealtvx = 2 * (precisionaltvx*recallaltvx) ./ (precisionaltvx + recallaltvx);

stats.simvoxel ={['True negative:']  [num2str(truenegative)];
                 ['False positive:'] [num2str(falsepositive)];
                 ['False negative:'] [num2str(falsenegative)];
                 ['True positive:']  [num2str(truepositive)];
                 ['Precision:']      [num2str(precision)];
                 ['Recall:']         [num2str(recall)]; 
                 ['F1 score:'] [num2str(f1score)]
                 ['Simulated voxel:'] [num2str(dipolevoxIN)]};        


stats.altvoxel ={['True negative:']  [num2str(truenegativealtvx)];
                 ['False positive:'] [num2str(falsepositive)];
                 ['False negative:'] [num2str(falsenegativealtvx)];
                 ['True positive:']  [num2str(truepositivealtvx)];
                 ['Precision:']      [num2str(precisionaltvx)];
                 ['Recall:']         [num2str(recallaltvx)]; 
                 ['F1 score:'] [num2str(f1scorealtvx)]
                 ['Alternative voxel:'] [num2str(maxvox)]};        
           
end


