% function phaseConsistencyWrapper

clear all;

% subject list - original batch
% all batches
slist = [1022:1026,1028:1029,1032, 1034:1040, 1042:1043, 1045:1052,1055,...
    1056:1061,1064:1066,1067,1070,1071,1073,1074,1075,1076,1077,1078,1079,1080,1081,1083];

roi = 'erc_right_park';

design_phA_on = true;

dataDir = '/gpfs/data/mnassar/lyu21/GridPlanet/';

% make directories for between  and within state testing
saveDir = [dataDir 'conditions_files_singleruns/stimMod_dropErrors_withinstate/'];

numRuns = 4;

numFolds = 6;

regressorIdx = [2,4,6,8];

%% get beta coefficients from each subject-run
for i = 1:length(slist)
    if slist(i) == 1066 || slist(i) == 1080
        numRuns = 3;
    else
        numRuns = 4;
    end
    for run = 1:numRuns
        [b]=marsbar_roi_singleRuns(slist(i),run, 'stimMod_dropErrors_state',roi);
        sin_cos_betas(i,run,:)=b([2,3,5,6,8,9,11,12]);
    end
end

%% loop subjects and runs - make conditions files
grid_betas = nan(length(slist),4,length(regressorIdx));

for i = 1:length(slist)
    
%     calculate mean betas for state 0 and state 1
    betas_state_0 = squeeze(mean(sin_cos_betas(:,:,[1,2,5,6]),2));
    betas_state_1 = squeeze(mean(sin_cos_betas(:,:,[3,4,7,8]),2));
    
    state0_idx = [1,2,5,6];
    state1_idx = [3,4,7,8];
    
    if slist(i) == 1066 || slist(i) == 1080
        numRuns = 3;
    else
        numRuns = 4;
    end
    
    for run = 1:numRuns
                    
            phi0=[];
            phi1=[];
            
            runList = 1:numRuns;
            runList = runList(runList~=run);
                
            sinPhi0 = sin_cos_betas(i,:,state0_idx(1)); %get the  sin & cos - for trials of interest
            cosPhi0 = sin_cos_betas(i,:,state0_idx(2));
            sinPhi1 = sin_cos_betas(i,:,state1_idx(1)); %get the  sin & cos - for trials of interest
            cosPhi1 = sin_cos_betas(i,:,state1_idx(2)); 

            phi0 = mod(atan2(sinPhi0,cosPhi0),pi*2);
            phi1 = mod(atan2(sinPhi1,cosPhi1),pi*2);
                
            mean_phi = [circ_mean(phi0,[],2),circ_mean(phi1,[],2)]./6;
            mean_phi_within = circ_mean([phi1(runList);phi0(runList)],[],2)/6;
                       
            %                 set up conditions file         
            SPM_designmaker_phaseConsistency_state_phi(slist(i), run, mean_phi_within, roi, saveDir, numFolds);   

    end
end
        
  %% run through GLMs and estimate betas in ROI

for i = 1:length(slist)
    if slist(i) == 1066 || slist(i) == 1080
        numRuns = 3;
    else
        numRuns = 4;
    end
    for run = 1:numRuns
        
        if design_phA_on
            %                 make SPM design file...
            spm_design_firstlvl_gridplanet_phA_state(slist(i), run,  roi, 'stimMod_phA_withinstate', true,true, numFolds);
        end
        
        %                 estimate using marsbar in our ROI
        [b]=marsbar_roi_phA_1run_deg(slist(i),run,'stimMod_phA_withinstate', roi);
        
        %                 keep betas for states 0 and 1, included and
        %                 excluded trials
        grid_betas(i,run,:) = b(regressorIdx);
        
        
        
        
    end
end

% average betas over runs 
mean_b=squeeze(nanmean(grid_betas,2));
% average betas over states 0/1 - for included
gmean_b(:,1) = nanmean(mean_b(:,[1,2]),2);
% for exluded trials
gmean_b(:,2) = nanmean(mean_b(:,[3,4]),2);

