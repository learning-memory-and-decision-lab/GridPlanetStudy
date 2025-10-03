function [b] = singleTrialRSA(subj,figon)
% load the concatenated single trial/miniblock data and calculate distances
% between miniblock pairs

% set the data path
data_path = '/gpfs/data/mnassar/lyu21/GridPlanet/';

% data path
roi_path = [data_path 'singleTrialROIData/spatialMeanPos_visual/'];

% load the data
load([roi_path sprintf('sub%d-roidata.mat',subj)]);

% find any nans in the data and remove those voxels
nanidx = (sum(isnan(trialByVoxel))>0);
trialByVoxel = trialByVoxel(:,~nanidx);

%mean center on each run
if subj == 1066 || subj == 1080
    numBlocks=3; %note, this is scan run not task block
    numtrial = 252;
else
    numtrial = 336;
    numBlocks = 4;
end
j=0;

for i = 1:numBlocks
    indicesRun = [j*(numtrial/numBlocks)+1:numtrial/numBlocks*i];
    trialByVoxel(indicesRun,:)=trialByVoxel(indicesRun,:)-mean(trialByVoxel(indicesRun,:));
    j=j+1;
end
subData=struct;
if subj == 1066 || subj == 1080
    subData=load([data_path sprintf('behavioralData/gridplanet_Scanning_subj_%d_block_3.mat',subj)]);
else
    subData=load([data_path sprintf('behavioralData/gridplanet_Scanning_subj_%d_block_4.mat',subj)]);
end

[T] = gridAnalysisFunction(subj,[1066,1080]);


% get regressors
[R, exclude] = makeHypothesisRDM_gridCircularDistance_v3(T,subData,[],0,0);


trialByVoxel = trialByVoxel(exclude==0,:);

% make a correlation distance matrix that is mb x mb
rdm = 1-corrcoef(trialByVoxel');
d = pdist(trialByVoxel,'correlation');
d = squareform(d);
if figon == true
    figure;imagesc(rdm);
end


% make matrix of regressors
regressors = [ones(length(R.runRegressor),1), ...
    R.runRegressor, R.blockRegressor, R.stateRegressor, ...
    R.movementShieldRegressor,...
    R.spatialMeanPositionRegressor, R.cogMeanPositionRegressor,...
    R.colorLABRegressor,...
    R.colorRegressor, R.colorWithinStateRegressor,...
    R.startShieldRegressor, R.endShieldRegressor,...
    R.transferTrialRegressor,R.angleChangeRegressor,...
    R.lagRegressors(:,:,1),...
    R.lagRegressors(:,:,2),R.lagRegressors(:,:,3),R.lagRegressors(:,:,4),...
    R.lagRegressors(:,:,5),R.lagRegressors(:,:,6),R.lagRegressors(:,:,7),...
    R.lagRegressors(:,:,8),R.lagRegressors(:,:,9),R.lagRegressors(:,:,10),...
    R.lagRegressors(:,:,11),R.lagRegressors(:,:,12),R.lagRegressors(:,:,13),...
    R.lagRegressors(:,:,14),R.lagRegressors(:,:,15),R.lagRegressors(:,:,16)];

% specify regressor names
regname = {'intercept','runs','blocks', 'state' ...
    'movementShield', ...
    'spatialMeanPos','cogMeanPos','colorLAB',...
    'color','colorWithinState', ...
    'startShield','endShield',...
    'transferTrial','angleChange',...
    'lag1','lag2', 'lag3', 'lag4', ...
    'lag5', 'lag6', 'lag7', 'lag8', 'lag9', 'lag10', 'lag11', 'lag12', 'lag13', ...
    'lag14', 'lag15', 'lag16'};


% get lower triangle
mask=tril(true(size(d)),-1);
dVec=d(mask);

% regress predictors onto the data!
b = regress(dVec,regressors);
disp('jee')