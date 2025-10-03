% this script runs the between-state phase consistency analysis for a given
% ROI

clear all;

% subject list
% all batches
slist = [1022:1026,1028:1029,1032, 1034:1040, 1042:1043, 1045:1052,1055,...
    1056:1061,1064:1066,1067,1070,1071,1073,1074,1075,1076,1077,1078,1079,1080,1081,1083];

% specify roi here
roi = '';

design_phA_on = true;

dataDir = '/gpfs/data/mnassar/lyu21/GridPlanet/';

% make directories for between  and within state testing
saveDir = [dataDir 'conditions_files_allruns/stimMod_dropErrors_state/'];

% index for regressors to use
regressorIdx = [2,4,6,8];

%% get beta coefficients from each subject-run
for i = 1:length(slist)
    [b]=marsbar_roi_allRuns(slist(i), 'stimMod_dropErrors_state',roi);
    sin_cos_betas(i,:)=b([2,3,5,6,8,9,11,12]);
end
%% load beta coeff if prev step is already complete
for i = 1:length(slist)
    load(fullfile(dataDir,'allrunsGLMs/stimMod_dropErrors_state',sprintf('%d/betas_dropErrors_%s.mat',slist(i),roi)));
    sin_cos_betas(i,:)=b([2,3,5,6,8,9,11,12]);
end
%% loop subjects and runs - make conditions files
for i = 1:length(slist)
    
    phi0=[];
    phi1=[];

    numFolds = 6;

    phi0 = mod(atan2(sin_cos_betas(i,1),sin_cos_betas(i,2)),pi*2)./numFolds;
    phi1 = mod(atan2(sin_cos_betas(i,3),sin_cos_betas(i,4)),pi*2)./numFolds;


    %                 set up conditions file        
    SPM_designmaker_phaseConsistency_state_allruns_phi(slist(i), [phi0,phi1],roi, saveDir, numFolds);   

end
        
  %% run through GLMs and estimate betas in ROI
grid_betas = zeros(length(slist),length(regressorIdx));
for i = 1:length(slist)

    %                 make SPM design file...
    spm_design_firstlvl_gridplanet_phA_state_allruns(slist(i),  roi, 'stimMod_phA_state', true, numFolds,0);

    %                 estimate using marsbar in our ROI
    [b]=marsbar_roi_phA_allruns_deg(slist(i),'stimMod_phA_state', roi);

    %                 keep betas for states 0 and 1, included and
    %                 excluded trials
    grid_betas(i,:) = b(regressorIdx);

    
end


% average betas over states 0/1 - for included
gmean_b(:,1) = nanmean(grid_betas(:,[1,2]),2);
% for excluded trials
gmean_b(:,2) = nanmean(grid_betas(:,[3,4]),2);



