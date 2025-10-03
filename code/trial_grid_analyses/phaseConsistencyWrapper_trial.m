% function phaseConsistencyWrapper

clear all;



% all batches
slist = [1022:1026,1028:1029,1032, 1034:1040, 1042:1043, 1045:1052,1055,...
    1056:1061,1064:1066,1067,1070,1071,1073,1074,1075,1076,1077,1078,1079,1080,1081,1083];


roi = 'ppc';

dataDir = '/gpfs/data/mnassar/lyu21/GridPlanet/';

% specify where SPM mat files are for loading data and regressors
spmMatDir = fullfile(dataDir,'singlerunsGLMs','stimMod_dropErrors_state');

numRuns = 4;

numFolds = 6;

% specify TR in sec
TR = 2;

regressorIdx = [2,4,6,8];

%% get beta coefficients from each subject-run
for i = 1:length(slist)
    if slist(i) == 1066 || slist(i) == 1080
        numRuns = 3;
    else
        numRuns = 4;
    end
    for run = 1:numRuns
        [b,roi_file]=marsbar_roi_singleRuns(slist(i),run, 'stimMod_dropErrors_state',roi);
        sin_cos_betas(i,run,:)=b([2,3,5,6,8,9,11,12]);
    end
end

%% loop subjects and runs - make conditions files
grid_betas = nan(length(slist),4,length(regressorIdx));
mean_phi0 = nan(length(slist),4);
mean_phi1 = nan(length(slist),4);

avgPhiOn = 1;
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

            mean_phi0(i,run) = circ_mean(phi0(runList),[],2)./6;
            mean_phi1(i,run) = circ_mean(phi1(runList),[],2)./6;
                                       
    end
end

%% trial classification

% make a cell array to keep classifcation data
state_class = struct;
state_acc = struct;

% make nan arrays to keep accuracy and classification info for each subj x
% run
state_class_m = nan(length(slist),4);
state_class_acc_m = nan(length(slist),4);
state_class_acc_error_m = nan(length(slist),4);
state_class_acc_corr_m = nan(length(slist),4);
state_class_acc_corrAlt_m = nan(length(slist),4);

% loop through participants
for i  = 1:length(slist)

    
    if slist(i) == 1066 || slist(i) == 1080
        numRuns = 3;
    else
        numRuns = 4;
    end


        % load raw behavioral data
    if slist(i) == 1066 || slist(i) == 1080
        subData=load([dataDir sprintf('behavioralData/gridplanet_Scanning_subj_%d_block_3.mat',slist(i))]);
    else
        subData=load([dataDir sprintf('behavioralData/gridplanet_Scanning_subj_%d_block_4.mat',slist(i))]);
    
    end

    % get subject data info
    [T, store] = gridAnalysisFunction(slist(i),[1066,1080]);

    XCenter = subData.xCenter;
    YCenter = subData.yCenter;
  
    %get error calculated from center of circle
    T.idealAngle = atan2d((T.MeanY-YCenter),(T.MeanX-XCenter));
    T.actualAngle = atan2d((T.RespY-YCenter),(T.RespX-XCenter));
    T.idealAngle  = mod(T.idealAngle+360,360);
    T.actualAngle  = mod(T.actualAngle+360,360);
    T.angleDiff = abs(rad2deg(circ_dist(deg2rad(T.idealAngle),deg2rad(T.actualAngle))));

    %      calculate distance travelled 
    T.travelDist = sqrt(((T.ShieldStartX-T.RespX).^2)+(T.ShieldStartY-T.RespY).^2);

     %     calculate error for alternative state
    T.idealAngleAlt = zeros(length(T.idealAngle),1);
    for k = 1:height(T)
         if mod(T.Block(k),2) == 0
             T.idealAngleAlt(k,1) = mod(T.idealAngle(k) - 90,360);
         else
             T.idealAngleAlt(k,1) = mod(T.idealAngle(k) + 90,360);
         end
    end

    T.angleDiffAlt = abs(rad2deg(circ_dist(deg2rad(T.idealAngleAlt),deg2rad(T.actualAngle))));
     
%    trials where angular distance is not in line with either state
    angleThresh = 45;
    moveThresh = 10;
    T.noMovement = T.travelDist <= moveThresh;
    T.errorExclude = T.angleDiff>angleThresh & T.angleDiffAlt > angleThresh & T.noMovement == 0;  
    T.errorExclude = T.angleDiff>angleThresh & T.angleDiffAlt > angleThresh & T.noMovement == 0;
    T.correctAngle = T.angleDiff<angleThresh & T.angleDiffAlt >=angleThresh & T.noMovement == 0;
    T.correctAngleAlt = T.angleDiff>=angleThresh & T.angleDiffAlt<angleThresh & T.noMovement == 0;

    % initialize rho coefficient matrix
    rho0.subj{i} = nan(numRuns,height(T)/numRuns);
    rho1.subj{i} = nan(numRuns,height(T)/numRuns);

    % loop through runs 
    for run = 1:numRuns
        

        % % first, we are going to load in the filtered and whitened timeseries 
        % data for our ROI
        spmSubjPath = fullfile(spmMatDir,mat2str(slist(i)),mat2str(run),'SPM.mat');
        [y,design_obj] = timeseries_extraction(spmSubjPath,roi_file);

        % select run data 
        runT = T(T.MRIBlock==run,:);

%       make matrix to keep data
        state_class(i).run{run} = nan(height(runT),1);
        state_acc(i).run{run} = nan(height(runT),1);

        % get event timings for trials in run
        runTrialTimes = store(run).Time(contains(store(run).Event,'sim'))-store(run).Time(strcmp(store(run).Event,'BLOCKSTART'));

        % load the SPM mat file
        load(spmSubjPath);

        % get regressors from SPM mat - removing existing parametric
        % modulators
        idx_nuis_reg = ~[contains(SPM.xX.name,'bf') | contains(SPM.xX.name,'Sin') | contains(SPM.xX.name,'Cos')];
        X = SPM.xX.xKXs.X(:,idx_nuis_reg);
        idx_reg = contains(SPM.xX.name,'state') & ~[contains(SPM.xX.name,'Sin') | contains(SPM.xX.name,'Cos')];
        X = [X,sum(SPM.xX.xKXs.X(:,idx_reg),2)];

        % regress these out from data
        b_X = regress(y,X);
        y_hat = X*b_X;

        % calculate residuals
        y_res = y - y_hat;
        % get trial timings in TRs
        runTrialTimesTR = runTrialTimes/TR;

        %  estimate responses based on phi and travelled angle
        modState0 = cos(numFolds*(deg2rad(runT.respAngle)-mean_phi0(i,run)));
        modState1 = cos(numFolds*(deg2rad(runT.respAngle)-mean_phi1(i,run)));

        modState0 = modState0 - mean(modState0);
        modState1 = modState1 - mean(modState1);

        modState0_cell{i,run} = modState0;
        modState1_cell{i,run} = modState1;

        % now loop through trials...!
        for t = 1:height(runT)

            start_idx = runTrialTimesTR(t);

            % get vector for TRs in trial to add regressor depending on
            % duration
            end_idx = runTrialTimesTR(t)+(runT.WaitSimTime(t)/TR);

            % make convolved BOLD response
            trial_bold = sim_bold(TR,start_idx,end_idx,1,length(y));
            trial_mod0 = sim_bold(TR,start_idx,end_idx,modState0(t),length(y));
            trial_mod1 = sim_bold(TR,start_idx,end_idx,modState1(t),length(y));

            % apply design filter to both modulators
            trial_mod0_marsy = marsy(trial_mod0); 
            trial_mod0_filt = apply_filter(design_obj,trial_mod0_marsy);
            trial_mod_0 = summary_data(trial_mod0_filt);

            trial_mod1_marsy = marsy(trial_mod1); 
            trial_mod1_filt = apply_filter(design_obj,trial_mod1_marsy);
            trial_mod_1 = summary_data(trial_mod1_filt);

            % % get correlation between regressors
            rho0.subj{i}(run,t) = corr(y_res,trial_mod0,'type','pearson'); 
            rho1.subj{i}(run,t) = corr(y_res,trial_mod1,'type','pearson'); 


            if rho0.subj{i}(run,t) ~= rho1.subj{i}(run,t) 
                if rho0.subj{i}(run,t) > rho1.subj{i}(run,t)
                    state_class(i).run{run}(t) = 0;
                else
                    state_class(i).run{run}(t) = 1;
                end
            else
                state_class(i).run{run}(t) = NaN;
            end
        end

        % calculate class accuracy
        stateBlock = 1-mod(runT.Block,2);
        state_acc(i).run{run} = double(state_class(i).run{run} == stateBlock);
        state_acc(i).run{run}(isnan(state_class(i).run{run})) = NaN;

        state_class_m(i,run) = mean(state_class(i).run{run},'omitnan');
        state_class_acc_m(i,run) = mean(state_acc(i).run{run},'omitnan');

        state_class_acc_error_m(i,run) = mean(state_acc(i).run{run}(runT.errorExclude == 1),'omitnan');

        state_class_acc_corr_m(i,run) = mean(state_acc(i).run{run}(runT.correctAngle == 1),'omitnan');
        state_class_acc_corrAlt_m(i,run) = mean(state_acc(i).run{run}(runT.correctAngleAlt == 1),'omitnan');

        % count non-nan trials
        state_class_nonnan(i,run) = sum(~isnan(state_class(i).run{run}));
    end
end



