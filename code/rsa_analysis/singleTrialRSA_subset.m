function [b] = singleTrialRSA_subset(subj,figon)
% load the concatenated single trial/miniblock data and calculate distances
% between miniblock pairs

% set the data path
data_path = '/gpfs/data/mnassar/lyu21/GridPlanet/';

% data path
roi_path = [data_path 'singleTrialROIData/spatialMeanPos_ofc_bi/'];

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


% load behavioral data
load(fullfile(data_path,'behavioralData',sprintf('extracted_data/%d_extracted_data_ses_2.mat',subj)));

% make hypothesis RDMS
[R, exclude] = makeHypothesisRDM_gridCircularDistance_v3(T,subData,[],0,0);


% trialByVoxel = trialByVoxel(exclude==0,:);

% choose trials to include in analysis
idx = T.transfer_trial == 1 & T.novel_color==0 & T.Block > 1;

% select subset of trials from trial x voxel matrix
trialByVoxel = trialByVoxel(idx,:);

% make a correlation distance matrix that is mb x mb
rdm = 1-corrcoef(trialByVoxel');
d = pdist(trialByVoxel,'correlation');
d = squareform(d);

if figon == true
    figure;imagesc(rdm);
end

% concatenate dregressors
regressors = [ones(length(R.runRegressor),1),...
    R.spatialMeanPositionRegressor];

% get lower triangle of distance matrix
mask=tril(true(size(d)),-1);
dVec=d(mask);

% filter out regressors
for j = 1:size(regressors,2)
    rdm_reg = squareform(regressors(:,j),'tomatrix');
    rdm_reg_reduced = rdm_reg(idx,idx);
    regressors_reduced(:,j) = squareform(rdm_reg_reduced,'tovector');
end

% regress predictors onto the data!
b = regress(dVec,regressors_reduced);
disp('jee')