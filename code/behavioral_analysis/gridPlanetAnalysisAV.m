function gridPlanetAnalysisAV(slist,session)
% this function runs main analyses of behavioral data for grid planet study
% requires a list of subject numbers (below) and session number to run.
% also produces main figures for the manuscript
% Yu et al. 2025

% sessions 1+2
% slist = [1022	1023	1024	1025	1026	1028	1029	1032	1034	1035	1036	1037	1038	1039	1040	1042	1043	1045	1046	1047	1048	1049	1050	1051	1052	1055	1056	1057	1058	1059	1060	1061	1064	1065	1066	1067	1070	1071	1073	1074	1075	1076	1077	1078	1079	1080	1081	1083];

% session 3
% slist = [1022	1023	1024	1025	1026	1028	1029	1032	1034	1035	1036	1037	1038	1039	1040	1042	1043	1046    1047	1048	1049	1050	1051	1052	1055	1056	1057	1058	1059	1060	1061	1064	1065	1066	1067	1070	1071	1073	1074	1075	1076	1077	1078	1079	1080	1081	1083];

% specify where data lives here
datadir = '';

% initialize transfer error var
transfer_error_block = nan(length(slist),8);
repeat_error_block = nan(length(slist),8); 
transferTrialsErrorM = nan(length(slist),4);
stateChangeTrials = nan(length(slist),8,5);
if session == 2
    nreps = 6;
    errorBlockRep = nan(length(slist),nreps,8);
    errorBlockRep_novel = nan(length(slist),nreps,8);
else
    nreps = 8;
    errorBlockRep = nan(length(slist),nreps,8);
    errorBlockRep_novel = nan(length(slist),nreps,8);
end


% run through list of subjects and load data for appropriate session
for i = 1:length(slist)
    
    %     load the center of the screen info from the raw data (yes, this is
    %     stupid)
    if session == 1 
        load(fullfile(datadir,sprintf('Data/session%d/rawData/gridplanet_V2_subj_%d.mat',session,slist(i))),'xCenter','yCenter');
    elseif session == 2
        if slist(i) == 1066 || slist(i) == 1080
            load(fullfile(datadir,sprintf('Data/fMRIdata/gridplanet_Scanning_subj_%d_block_3.mat',slist(i))),'xCenter','yCenter');
        else
            load(fullfile(datadir,sprintf('Data/fMRIdata/gridplanet_Scanning_subj_%d_block_4.mat',slist(i))),'xCenter','yCenter');
        end
    elseif session == 3
        load(fullfile(datadir,sprintf('Data/session%d/rawData/gridplanet_V2_subj_%d_angleswitch.mat',session,slist(i))),'xCenter','yCenter');
    end
    %     load the data
    load(fullfile(datadir,sprintf('Data/session%d/sortedData_subject%d_session%d.mat',session,slist(i),session)));
    
    %get error
    distMean =  sqrt((T.MeanX-T.RespX).^2 + (T.MeanY - T.RespY).^2); 
    
    %     calculate correct angle for each trial
    %     calculated from center of circle
    idealAngle = atan2d((T.MeanY-yCenter),(T.MeanX-xCenter));
    actualAngle = atan2d((T.RespY-yCenter),(T.RespX-xCenter));
    idealAngle  = mod(idealAngle+360,360);
    actualAngle  = mod(actualAngle+360,360);
    T.angleDiff = abs(rad2deg(circ_dist(deg2rad(idealAngle),deg2rad(actualAngle))));
    T.distMeanEuc = distMean;
    T.actualAngle = actualAngle;
    T.idealAngle = idealAngle;
    
%     specify MRI runs
    T.MRIblock(T.Block==1 | T.Block==2) = 1;
    T.MRIblock(T.Block==3 | T.Block==4) = 2;
    T.MRIblock(T.Block==5 | T.Block==6) = 3;
    T.MRIblock(T.Block==7 | T.Block==8) = 4;
    

    if session == 2
        lenScanBlk = 84;
    elseif session == 1 || session == 3
        lenScanBlk = 80;
    end
    
    %for dropping trials analysis:
    state2 = zeros(lenScanBlk,1);
    state2(lenScanBlk/2+1:end) = 1;
    if session == 2 && (slist(i) == 1066 || slist(i) == 1080)
        state=repmat(state2,3,1);
    else
        state=repmat(state2,4,1);
    end
    T.state = state;


    % record state for first and second presentation of each novel color so
    % that we can know its initial state and transfer state
    NovelColorStateOrder = nan(13,2);
    % make a vector for noting transition trials
    transition_trial=zeros(height(T),1);
    for j = 1:height(T)
    %     identify novel colors and save their initial and transfer states...
        if T.ColorID(j) > 5
            if isnan(NovelColorStateOrder(T.ColorID(j),1))
                NovelColorStateOrder(T.ColorID(j),1)=T.state(j);
            elseif T.state(j) ~=  NovelColorStateOrder(T.ColorID(j),1)
                 NovelColorStateOrder(T.ColorID(j),2) = T.state(j);
            end
        end
        % identify transfer trials
        if j > 1
            if T.state(j) ~= T.state(j-1)
                transition_trial(j) = 1;
            end
        else
            transition_trial(j) = 1;
        end
    end

    % label those transfer trials
    transition_idx = find(transition_trial);
    state_change_trial = zeros(height(T),1);
    transfer_trial = zeros(height(T),1);
    for j = 1:length(transition_idx)
    %     label transfer trials for old colors
        transfer_trial(transition_idx(j)+1:transition_idx(j)+4)=1;
        if session == 2
            state_change_trial(transition_idx(j):transition_idx(j)+6)=1;
        else
            state_change_trial(transition_idx(j):transition_idx(j)+4)=1;
        end
    end

    %     determine whether novel color trial is new or transfer trial
    for j = 1:height(T)
        if T.ColorID(j) > 5
            if T.state(j) == NovelColorStateOrder(T.ColorID(j),2) && state_change_trial(j) == 1
                transfer_trial(j) = 1;
            end
        end
    end

    % append this new info to table T
    T.transfer_trial = transfer_trial;
    T.state_change_trial = state_change_trial;
    T.novel_color = double(T.ColorID>5);
    T.novel_color_ltm = T.ColorID==13;


    states = unique(T.state);
    colors = unique(T.ColorID);


%   if this is session 1 or 2, the angular shift between states is always
%   90 degrees. However, if this is session 3, then the angular shift is
%   changing on each block and we need to calculate 
    if session == 1 || session == 2
        angleShift = 90;
    elseif session == 3
        blockNums = unique(T.Block);
        colorNums = unique(T.ColorID);
        
        for b = blockNums'
            for c = colorNums'
                colorMeanX_by_Block(c,b) = mean(T.MeanX(T.ColorID==c & T.Block==b));
                colorMeanY_by_Block(c,b) = mean(T.MeanY(T.ColorID==c & T.Block==b));
            end
        end
%         now run through the series of blocks and calculate angular
%         distance from center
        for b = 1:length(blockNums)
            angleBlock(b,1) = mod(atan2d((colorMeanY_by_Block(1,b)-yCenter),(colorMeanX_by_Block(1,b)-xCenter)),360);
        end
%         now calculate the angular distance between blocks
        for b = 2:length(blockNums)
            angleChangeBlock(b,1) = mod(rad2deg(circ_dist(deg2rad(angleBlock(b)),deg2rad(angleBlock(b-1)))),360);
        end
    end
                
            
        
    %     calculate error for alternative state
        idealAngleAlt = zeros(length(idealAngle),1);
        for k = 1:height(T)
%             if this is session 3, then flexibly update angleShift by
%             block number
            if session == 3
                if T.Block(k) == 1
                    angleShift = 90;
                else
                    angleShift = angleChangeBlock(T.Block(k));
                end
            end

            if session == 1 || session == 2
                 if T.Block(k) == 1 || T.Block(k) == 3 || T.Block(k) == 5 || T.Block(k) == 7 
                     idealAngleAlt(k,1) = mod(T.idealAngle(k) + angleShift,360);
                 else
                     idealAngleAlt(k,1) = mod(T.idealAngle(k) - angleShift,360);
                 end
            elseif session == 3
                idealAngleAlt(k,1) = mod(T.idealAngle(k) - angleShift,360);
            end
                
        end
    T.idealAngleAlt = idealAngleAlt;
    T.angleDiffAlt = abs(rad2deg(circ_dist(deg2rad(idealAngleAlt),deg2rad(actualAngle))));
    T.normError = T.angleDiff./(T.angleDiff+T.angleDiffAlt);
%     
   

% first, let's compare the norm error on new colors to the error on
% transfer trials to the normed error on the first transfer trial...

    if session == 2
        newColorTransferError(i,1) = mean(T.normError(T.novel_color==1 & T.novel_color_ltm==0 & T.transfer_trial==1 & T.Block ~=1));
    end
    firstTrialTransferError(i,1) = mean(T.normError(T.novel_color==0 & T.transfer_trial==0 & T.state_change_trial == 1 & T.Block ~=1));
    
    
%     get the average error for each transfer trial - in order
    transferTrialsError = T.normError(T.transfer_trial == 1 & T.novel_color == 0 & T.Block ~=1);
        
%     get the error for each color on each block
    k = 0;
    for t = 1:4
        transferTrialsErrorM(i,t) = mean(transferTrialsError([t:4:length(transferTrialsError)]));
    end
%     cycle through the 'blocks' and get the average transfer error for
%     each
    blockNums = unique(T.Block);
    stimReps = unique(T.StimulusRep);
    
    for b = blockNums'
        transfer_error_block(i,b) = mean(T.normError(T.Block == b & T.transfer_trial == 1 & T.state_change_trial == 1 & T.novel_color == 0));
        repeat_error_block(i,b) = mean(T.normError(T.Block == b & T.transfer_trial == 0 & T.state_change_trial == 0 & T.novel_color == 0));
    end
    
%   measure the error for each stimulus repetition at each block
    for b = blockNums'
        for s = stimReps'
            errorBlockRep(i,s,b) = mean(T.normError(T.Block == b & T.StimulusRep == s & T.novel_color == 0));
%         get seperate measures for novel colors in session 2  
            if session == 2
                errorBlockRep_novel(i,s,b) = mean(T.normError(T.Block == b & T.StimulusRep == s & T.novel_color == 1 & T.novel_color_ltm==0));
            end
        end
    
%     get the error for the first set of state change trials at each block
       stateChangeTrials(i,b,:) = T.normError(T.Block == b & T.state_change_trial == 1 & T.novel_color == 0);

    end
    
%     fill out the trial number since each state change
    blocksAll = unique(T.Block);
    for k = 1:length(blocksAll)
        T.trialFromStateChange(T.Block == blocksAll(k)) = 1:sum(T.Block==blocksAll(k));
    end
    
    
    
    
% % % % --- Statistical analysis ----- % % % % %  %
    idx = T.novel_color == 0 & T.Block ~=1;
    [beta_repreg(i,:)] = regress(T.normError(idx),[ones(height(T(idx,:)),1),zscore(T.MRIblock(idx)),zscore(T.StimulusRep(idx)),zscore(zscore(T.MRIblock(idx)).*zscore(T.StimulusRep(idx)))]);

    idx = T.state_change_trial == 1 & T.novel_color == 0 & T.Block ~=1;
    [beta_trialreg(i,:)] = regress(T.normError(idx),[ones(height(T(idx,:)),1),zscore(T.MRIblock(idx)),zscore(T.trialFromStateChange(idx)),zscore(zscore(T.MRIblock(idx)).*zscore(T.trialFromStateChange(idx)))]);

    
    idx = T.state_change_trial == 1 & T.novel_color == 0 & T.Block ~=1 & T.trialFromStateChange == 2;
    [beta_trial_block_2(i,:)] = regress(T.normError(idx),[ones(height(T(idx,:)),1),zscore(T.MRIblock(idx))]);

    idx = T.state_change_trial == 1 & T.novel_color == 0 & T.Block ~=1 & T.trialFromStateChange == 5;
    [beta_trial_block_5(i,:)] = regress(T.normError(idx),[ones(height(T(idx,:)),1),zscore(T.MRIblock(idx))]);
    

end

%     ###################################################################
%   ok, now let's start doing our analyses.


% confidence interval
ci = 0.975;

% set font size
fsize = 16;

% marker size
markerSz = 80;

% Set color vector
cStart = [0.05, 0.3, 0.35];
cEnd =   [0.2, 0.8, 0.6];
color_transfer_trials = interp1([1;4],[cStart;cEnd],(1:4)');



% make a figure showing the change in the transfer error from
% block-to-block - excluding the first one
figure;
transfer_error_block_m = mean(transfer_error_block,'omitnan');
transfer_error_block_ci = tinv(ci,length(slist)-1).*(std(transfer_error_block,'omitnan')./sqrt(length(slist)));
shadedErrorBar([1:7],transfer_error_block_m(2:8),transfer_error_block_ci(2:8),{'color',[0.4,0.05,0.75]});
xticks(1:7);
xlim([1,7]);
ylim([0,1]);
xlabel('state change');
ylabel('Normalized Transfer Error');
set(gca,'fontsize',fsize);
set(gcf, 'Position',  [100, 100, 300, 400])

%     make figure for old color transfer vs first trial state change error
figure;
% for t = 1:4
scatter(firstTrialTransferError,mean(transferTrialsErrorM,2),markerSz,'filled','MarkerEdgeColor',[0.1 0.3 .6],'MarkerFaceColor',[0.2 0.5 .8], 'LineWidth',1.5); hold all;
% end
ylim([0,1]);
xlim([0,1]);
plot([0,1],[0,1],'k--','LineWidth',3);
box on;
xlabel('First Trial State Change Error');
ylabel('Transfer Error');
set(gca,'fontsize',fsize);
set(gcf, 'Position',  [100, 100, 325, 300])

% test if error on transfer trials is below error on first trial
transferErrorDiff = mean(transferTrialsErrorM,2) - firstTrialTransferError;
[~,p_transferdiff,~,stats_transferdiff] = ttest(transferErrorDiff);
[~,p_transferdiff_chance,~,stats_transferdiff_chance] = ttest(mean(transferTrialsErrorM,2),0.5);

[~,p_transferdiff_chance_b1,~,stats_transferdiff_chance_b1] = ttest(transferTrialsErrorM(:,1),0.5);
[~,p_transferdiff_chance_b4,~,stats_transferdiff_chance_b4] = ttest(transferTrialsErrorM(:,4),0.5);


% compare error on first trial to chance
[~,p_firstTransfer_chance,~,stats_firstTransfer_chance] = ttest(firstTrialTransferError,0.5);

    
% session 2 plots only
    if session == 2
    %     make figure for new color - first trial transfer error
%         figure;
        hold all;
        scatter(firstTrialTransferError,newColorTransferError,markerSz,'filled','MarkerEdgeColor',[0.5 0 .5],'MarkerFaceColor',[0.7 0 .7], 'LineWidth',1.5); hold all;
        ylim([0,1]);
        xlim([0,1]);
        plot([0,1],[0,1],'k--','LineWidth',3);
        box on;
        xlabel('First Trial State Change Error');
        ylabel('New Color Transfer Error');
        set(gca,'fontsize',fsize);
%         set(gcf, 'Position',  [100, 100, 450, 400])
        
        % test if error on transfer trials on novel colors is below error on first trial
        transferErrorDiffNovel = newColorTransferError - firstTrialTransferError;
        [~,p_transferdiff_novel,~,stats_transferdiff_novel] = ttest(transferErrorDiffNovel);
        [~,p_transferdiff_novel_chance,~,stats_transferdiff_novel_chance] = ttest(newColorTransferError,0.5);


    end
    
    
% make a figure for the change in the normalized transfer error across both
% blocks and stimulus repetitions


%  interpolate colors across 8 blocks...
cstart = [0.2,0.15,0.4];
cend = [0.85,0.3,0.65];
for c = 1:3
    cspace(:,c) = linspace(cstart(c),cend(c),8);
end

figure;
statechange_b = {[2];[3,4];[5,6];[7:8]};
for b = 1:size(statechange_b,1)
    blockrep_error_m = mean(mean(errorBlockRep(:,:,statechange_b{b}),3,'omitnan'),'omitnan');
    n_sb = sum(~isnan(errorBlockRep(:,:,b)));
    t = tinv(ci,n_sb);
    blockrep_error_ci = t.*(std(mean(errorBlockRep(:,:,statechange_b{b}),3,'omitnan'),'omitnan')./sqrt(n_sb));
    shadedErrorBar([1:nreps],blockrep_error_m,blockrep_error_ci,{'color',cspace(b,:)});       
    hold all;
end

ylim([0,1]);
xlim([1,length(stimReps)]);
xticks(stimReps);
xlabel('Stimulus Repetition');
ylabel('Normalized Error');
set(gca,'fontsize',fsize);
set(gcf, 'Position',  [100, 100, 450, 400])

%  ---- make a plot of the mean norm error across state reps over all
%  blocks
figure;
error_m = mean(mean(errorBlockRep,3,'omitnan'),'omitnan');
n_sb = length(slist);
t = tinv(ci,n_sb);
error_ci = t.*(std(mean(errorBlockRep,3,'omitnan'),'omitnan')./sqrt(n_sb));
shadedErrorBar([1:nreps],error_m,error_ci,{'color',cspace(b,:)});       
yyaxis right
xlim([1,length(stimReps)]);
xticks(stimReps);
xlabel('Stimulus Repetition');
ylabel('Normalized Error');
set(gca,'fontsize',fsize);
set(gcf, 'Position',  [100, 100, 450, 400])
yyaxis left
yticks([]);
ylim([0,1]);


%  ------ do the same for the novel trials -----
if session == 2
    figure;
    statechange_b = {[2];[3,4];[5,6];[7:8]};
    for b = 1:size(statechange_b,1)
        blockrep_error_m = mean(mean(errorBlockRep_novel(:,:,statechange_b{b}),3),'omitnan');
        n_sb = sum(~isnan(errorBlockRep(:,:,b)));
        t = tinv(ci,n_sb-1);
        blockrep_error_ci = t.*(std(mean(errorBlockRep_novel(:,:,statechange_b{b}),3),'omitnan')./sqrt(n_sb));
        shadedErrorBar([1:nreps],blockrep_error_m,blockrep_error_ci,{'color',cspace(b,:)});       
        hold all;
    end

    ylim([0,1]);
    xlim([1,length(stimReps)]);
    xticks(stimReps);
    xlabel('Stimulus Repetition');
    ylabel('New Color Normalized Error');
    set(gca,'fontsize',fsize);
    set(gcf, 'Position',  [100, 100, 450, 400])
end


% make figure for state change trials ------
figure;
statechange_b = {[2];[3,4];[5,6];[7:8]};
for b = 1:size(statechange_b,1)
   stateChangeError_m(b,:) = squeeze(mean(mean(stateChangeTrials(:,statechange_b{b},:),2,'omitnan'),'omitnan'));
   n_sb = squeeze(sum(~isnan(stateChangeTrials(:,b,:))));
   t = tinv(ci,n_sb-1);
   stateChangeError_ci(b,:) = t.*(squeeze(std(mean(stateChangeTrials(:,statechange_b{b},:),2,'omitnan'),'omitnan'))./sqrt(n_sb));

    % make a figure
    shadedErrorBar([1:5],stateChangeError_m(b,:),stateChangeError_ci(b,:),{'color',cspace(b,:)});       
    hold all;
end
xticks(1:5);
ylim([0,1]);
xlim([1,5]);
xlabel('Trials');
ylabel('Normalized Error')
set(gca,'fontsize',fsize);
set(gcf, 'Position',  [100, 100, 450, 400])

% make a figure for mean change across blocks for error after a state
% change
figure;
gmstateChangeError_m = squeeze(mean(mean(stateChangeTrials,2,'omitnan'),'omitnan'));
n_sb = length(slist);
t = tinv(ci,n_sb-1);
gmstateChangeError_ci = t.*(squeeze(std(mean(stateChangeTrials,2,'omitnan'),'omitnan'))./sqrt(n_sb));

% make a figure
shadedErrorBar([1:5],gmstateChangeError_m, gmstateChangeError_ci,{'color',cspace(b,:)});       
hold all;

xticks(1:5);
ylim([0,1]);
xlim([1,5]);
yyaxis right;
xlabel('Trials');
ylabel('Normalized Error')
set(gca,'fontsize',fsize);
set(gcf, 'Position',  [100, 100, 450, 400])


%     calculate residuals of logit transformed error measures
logitRE=log(mstaean(repeat_error_block,2)./(1-mean(repeat_error_block,2)));
logitTE=log(mean(transfer_error_block,2)./(1-mean(transfer_error_block,2)));
b = regress(logitTE,[ones(length(slist),1),logitRE]);
predTE =[ones(length(slist),1),logitRE]*b;
resTE=logitTE-predTE;


%  statistical analysis
[~,p_rep,~,stats_rep] = ttest(beta_repreg);
[~,p_trial,~,stats_trial] = ttest(beta_trialreg);
[~,p_trial_5,~,stats_trial_5] = ttest(beta_trial_block_5);
[~,p_trial_2,~,stats_trial_2] = ttest(beta_trial_block_2);






disp('jee');