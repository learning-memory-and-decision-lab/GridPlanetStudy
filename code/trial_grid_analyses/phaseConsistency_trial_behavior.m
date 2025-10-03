function phaseConsistency_trial_behavior
% this function tests the relationship between pearson correlation
% coefficients from single-trial modulators and behavioral measures of
% error in the grid planet task

% specify data directory
datadir = '/oscar/data/mnassar/lyu21/GridPlanet/';

behavioraldir = fullfile(datadir,'behavioralData/extracted_data');


% options for thresholding angle and movement for correct responses
    angle_thresh = 45;
    move_thresh = 10;

% boolean to decide which set of state data to use
withinstate = 0;

states = [0:1];

roi = 'pcc';


    if strcmp(roi,'ppc')
        load(fullfile(datadir,'roi_analysis/phaseConsistency_trial_betweenstate/ppc_trial_betweenstate.mat'));
    elseif strcmp(roi,'pcc')
        load(fullfile(datadir,'roi_analysis/phaseConsistency_trial_betweenstate/pcc_trial_betweenstate.mat'));
    elseif strcmp(roi,'rerc')
        load(fullfile(datadir,'roi_analysis/phaseConsistency_trial_betweenstate/rerc_trial_betweenstate.mat'));
    elseif strcmp(roi,'mpfc')
        load(fullfile(datadir,'roi_analysis/phaseConsistency_trial_betweenstate/mpfc_trial_betweenstate.mat'));
    elseif strcmp(roi,'rtpj')
        load(fullfile(datadir,'roi_analysis/phaseConsistency_trial_betweenstate/rtpj_trial_betweenstate.mat'));
    elseif strcmp(roi,'ltpj')
        load(fullfile(datadir,'roi_analysis/phaseConsistency_trial_betweenstate/ltpj_trial_betweenstate.mat'));
    elseif strcmp(roi,'lerc')
        load(fullfile(datadir,'roi_analysis/phaseConsistency_trial_betweenstate/lerc_trial_betweenstate.mat'));
    elseif strcmp(roi,'meanroi')
        load(fullfile(datadir,'roi_analysis/phaseConsistency_trial_betweenstate/meanroi_trial_betweenstate.mat'));
    end


% load coefficients from between state cross-validation
ppc_btstate = load('/oscar/data/mnassar/lyu21/GridPlanet/roi_analysis/phaseConsistency_state_allruns/ppc_gm_coggrid.mat');



% initialize mean coefficients for block-rep
rho_block_rep = nan(length(slist),4,6);

rho_full_all = nan(length(slist),336);

% load and process behavioral data
% loop through participants
for i  = 1:length(slist)

    % load behavioral data table
    load(fullfile(behavioraldir,sprintf('%d_extracted_data_ses_2.mat',slist(i))));



    
    %      calculate distance travelled 
    T.travelDist = sqrt(((T.ShieldStartX-T.RespX).^2)+(T.ShieldStartY-T.RespY).^2);

    if slist(i) == 1066 || slist(i) == 1080
        runNums = 3;
    else
        runNums = 4;
    end

    % load in the raw beahvioral data from scanner session...
    raw_data = load(fullfile(datadir,'behavioralData',sprintf('gridplanet_Scanning_subj_%d_block_%d.mat',slist(i),runNums)));
    
    % calculate curvature?
    max_curv = [];
    int_curv = [];
    initiation=[];
    direction360 = mod(raw_data.direction,360);
    for t = 1:height(T)
      max_curv(t,1) = max(abs(circ_dist(deg2rad(direction360(:,t)),deg2rad(T.respAngle(t)))))/T.travelDist(t);
      idx = direction360(:,t)~=0;
      int_curv(t,1) = mean(abs(circ_dist(deg2rad(direction360(idx,t)),deg2rad(T.respAngle(t)))))/T.travelDist(t);
      if sum(direction360(:,t))==0
          initiation(t,1) = NaN;
      else
          initiation(t,1) = find(direction360(:,t),1);
      end
    end
    
    m_max_curv(i,1) = mean(max_curv,'omitnan');
    m_int_curv(i,1) = mean(int_curv,'omitnan');

    rho_full = [];
    mod_full = [];

    % keep track of coefficients across runs
    rho_runs = [];
    rho_block = [];

    % concatenate rho from current state across runs
    for r = 1:runNums
        runT = T(T.MRIBlock==r,:);
        state_idx = mod(runT.Block,2)==0;
        if withinstate
            rho_run = [rho0.subj{i}(r,state_idx==0)';rho1.subj{i}(r,state_idx==1)'];
            rho_block(r,:) = mean([rho0.subj{i}(r,state_idx==0);rho1.subj{i}(r,state_idx==1)]);
        else
            rho_run = [rho1.subj{i}(r,state_idx==0)';rho0.subj{i}(r,state_idx==1)'];
            rho_block(r,:) = mean([rho1.subj{i}(r,state_idx==0);rho0.subj{i}(r,state_idx==1)]);
        end            
        rho_full = [rho_full;rho_run];
        rho_runs(r,:)=rho_run';

        % get modulators for each run and append together into a long
        % vector
        mod0_run = modState0_cell{i,r};
        mod1_run = modState1_cell{i,r};

        if withinstate
            mod_run = [mod0_run(state_idx==0);mod1_run(state_idx==1)];
        else
            mod_run = [mod1_run(state_idx==0);mod0_run(state_idx==1)];
        end            
        mod_full = [mod_full;mod_run];
    end

    %     fill out the trial number since each state change
    blocksAll = unique(T.Block);
    for k = 1:length(blocksAll)
        T.trialFromStateChange(T.Block == blocksAll(k)) = 1:sum(T.Block==blocksAll(k));
    end

    rho_subj(i,:) = mean(rho_runs);
    rho_subj_block(i,:) = mean(rho_block);

    % norm error transformation
    logitnorm = log(T.normError./(1-T.normError));
    rho_full = fisher_transform(rho_full);
    % rho_full = zscore(fisher_transform(rho_full));

    % calculate weights based on modulator absolute values
    mod_w = diag(abs(mod_full)/max(abs(mod_full)));
    mod_w2 = 1-mod_w;

    % get average of fisher pearson correlation coefficients across stim
    % repeats - for old colors
    T.noMovement = T.travelDist < move_thresh;
    for rep = 1:6
        rho_stimOld_rep(i,rep) = mean(rho_full(T.StimulusRep==rep & T.novel_color==0 & T.Block ~= 1));
        rho_stimNov_rep(i,rep) = mean(rho_full(T.StimulusRep==rep & T.novel_color==1 & T.Block ~= 1));
    end

    % average over the first few transfer trials
    for t = 1:5
        rho_transfer_trials(i,t) = mean(rho_full(T.trialFromStateChange ==t & T.novel_color==0 & T.Block ~= 1));
    end

    rho_corrOld_rep(i,1) = mean(rho_full(T.StimulusRep<=3 & T.novel_color==0 & T.angleDiff<angle_thresh & T.noMovement == 0));
    rho_corrOld_rep(i,2) = mean(rho_full(T.StimulusRep>3 & T.novel_color==0 & T.angleDiff<angle_thresh & T.noMovement == 0));

    rho_corrAltOld_rep(i,1) = mean(rho_full(T.StimulusRep<=3 & T.novel_color==0 & T.angleDiffAlt<angle_thresh & T.noMovement == 0));
    rho_corrAltOld_rep(i,2) = mean(rho_full(T.StimulusRep>3 & T.novel_color==0 & T.angleDiffAlt<angle_thresh & T.noMovement == 0));

    rho_errorOld_rep(i,1) = mean(rho_full(T.StimulusRep<=3 & T.novel_color==0 & T.angleDiffAlt >= angle_thresh & T.angleDiff>=angle_thresh & T.noMovement == 0));
    rho_errorOld_rep(i,2) = mean(rho_full(T.StimulusRep>3 & T.novel_color==0 & T.angleDiffAlt >= angle_thresh & T.angleDiff>=angle_thresh & T.noMovement == 0));

    % total trial type
    total_trial(i,1) = sum(T.angleDiff < angle_thresh & T.noMovement == 0);
    total_trial(i,2) = sum(T.angleDiffAlt < angle_thresh & T.noMovement == 0);
    total_trial(i,3) = sum(T.angleDiff >= angle_thresh & T.angleDiffAlt >= angle_thresh & T.noMovement == 0);
    total_trial(i,4) = sum(T.noMovement == 1);

    % get averages of fisher correlation coefficients for each repeat-state
    % change
    for b = 1:length(unique(T.MRIBlock))
        for rep = 1:6
            rho_block_rep(i,b,rep) = mean(rho_full(T.Block==b & T.StimulusRep==rep & T.novel_color==0));
        end
    end

    % get state numbers
    state_idx = mod(T.Block,2)==0;

    % get average for sliding window of angular error
    window = 30:30:180;
    windowSigned = -180:45:180;
    angleErrorAbs = (T.angleDiff);
    % calculate signed error
    angleErrorSigned = rad2deg(circ_dist(deg2rad(T.actualAngle),deg2rad(T.idealAngle)));
    % re-reference signed error so that rotations are consistent between states
    angleErrorRef = angleErrorSigned;
    angleErrorRef(state_idx==1) = angleErrorRef(state_idx==1)*-1;
        
    for w = 1:length(window)
        if w == 1
            aidx = angleErrorAbs<=window(w);
        else
            aidx = angleErrorAbs>window(w-1) & angleErrorAbs<=window(w);
        end
        rho_window(i,w) = mean(rho_full(aidx));
    end

    for w = 1:length(windowSigned)
        for s = 1:2
            if w == 1
                aidx = angleErrorSigned<=windowSigned(w) & state_idx == states(s);
            else
                aidx = angleErrorSigned>windowSigned(w-1) & angleErrorSigned<=windowSigned(w) & state_idx == states(s);
            end
            rho_window_signed(i,s,w) = mean(rho_full(aidx));
        end 
        if w == 1
            aidxRef = angleErrorRef<=windowSigned(w);
        else
            aidxRef = angleErrorRef>windowSigned(w-1) & angleErrorRef<=windowSigned(w);
        end
        rho_window_ref(i,w) = mean(rho_full(aidxRef));
    end

    % now, let's regress rho against behavior
    time_block = repmat([1:42,1:42]',runNums,1);
    % X = [ones(height(T),1),rho_full,T.WaitSimTime,T.StimulusRep,zscore(rho_full).*zscore(T.WaitSimTime)];
    X = [ones(height(T),1),rho_full];

    [b_rho(i,:)] = regress(logitnorm,X);
    [b_rho_early(i,:)] = regress(logitnorm(T.StimulusRep==1 | T.StimulusRep==2 | T.StimulusRep==3,1),X(T.StimulusRep==1 | T.StimulusRep==2 | T.StimulusRep==3,:));
    [b_rho_late(i,:)] = regress(logitnorm(T.StimulusRep==4 | T.StimulusRep==5 | T.StimulusRep==6,1),X(T.StimulusRep==4 | T.StimulusRep==5 | T.StimulusRep==6,:));
    Xtime = [X,(T.StimulusRep)];
    % Xtime = [Xtime,Xtime(:,2).*Xtime(:,3)];
    [b_rho_time(i,:)] = regress(logitnorm,Xtime);

    stimRepX = [ones(height(T),1),T.StimulusRep];
    stimrepidx = T.Block>=2;
    trialidx = [T.trialFromStateChange<=5 & T.Block>=2];
    trialX = [ones(height(T),1), T.trialFromStateChange];


    [b_stimrep(i,:)] = regress(rho_full(stimrepidx),stimRepX(stimrepidx,:));
    [b_trial(i,:)] = regress(rho_full(trialidx),trialX(trialidx,:));

    % test rho values on stim reps 
    for r = 1:5
        rho_mean(i,r) = mean(rho_full(T.StimulusRep == r & stimrepidx));
    end




    % regress against squared error
    X = [ones(height(T),1),rho_full,T.StimulusRep];
    [b_abs(i,:)] = regress(T.angleDiff,X);   

    %  behavioral regression against curviture of response
    Xcurv= [ones(height(T),1),T.StimulusRep,rho_full];

    b_curv(i,:) = regress(int_curv,Xcurv);

    b_rt(i,:) = regress(rho_full,[Xcurv(:,1),T.WaitSimTime]);

    % weighted regression
    bw_rho(i,:) = weighted_regress(logitnorm,X,mod_w);
    bw2_rho(i,:) = weighted_regress(logitnorm,X,mod_w2);

    % regression against pearson correlation
    b_time(i,:) = regress(rho_full,[ones(length(time_block),1),time_block]);
    b_rep(i,:) = regress(rho_full,[ones(height(T),1),T.StimulusRep,T.MRIBlock]);
    b_rep_beh(i,:) = regress(T.normError,[ones(height(T),1),T.StimulusRep]);

 
    
    rho_corr(i,1) = mean(rho_full(T.angleDiff<angle_thresh & T.noMovement == 0)); 
    rho_corrAlt(i,1) = mean(rho_full(T.angleDiffAlt<angle_thresh & T.noMovement == 0)); 
    rho_corrAltLate(i,1) = mean(rho_full(T.angleDiffAlt<angle_thresh & T.noMovement == 0 & T.StimulusRep > 2)); 
    rho_error_all(i,1) = mean(rho_full(T.angleDiff>angle_thresh & T.noMovement == 0)); 
    rho_error(i,1) = mean(rho_full(T.angleDiffAlt>angle_thresh & T.angleDiff>angle_thresh & T.noMovement == 0)); 
    rho_state_change(i,1) =  mean(rho_full(T.transfer_trial == 1  & T.noMovement == 0)); 

    rho_full_all(i,1:length(rho_full)) = rho_full;
end


plotson = true;
if plotson

    
    % set up colors for subjects
    cstart = [0.1,0.1,0.1];
    cend  = [0.5,0.2,0.9];
    for r = 1:3
        cspace(:,r) = linspace(cstart(r),cend(r),size(rho_block_rep,2));
    end


    % figure averaging across all blocks
    figure;
    m_rho_rep = mean(rho_stimOld_rep);
    sample_sz = size(rho_stimOld_rep,1);
    ci_rho_rep = [std(rho_stimOld_rep)./sqrt(sample_sz)] .* tinv(0.975,sample_sz-1); 
    errorbar(m_rho_rep,ci_rho_rep,'CapSize',0,'Color',repmat(0.8,1,3));
    hold all;plot(m_rho_rep,'b','LineWidth',2);
    xlim([0,7])
    xticks(1:6)
    xlabel('stim repeat');
    ylabel('mean Fisher correlation coeff')

    % make shaded error bar version of above
    figure;
    shadedErrorBar([1:6],m_rho_rep,ci_rho_rep);
    xlim([0,7])
    xticks(1:6)
    xlabel('stim repeat');
    ylabel('mean Fisher correlation coeff')

    % shaded error bar plot of rho across transfer trials for old colors
    figure;
    m_rho_transfer_trial = mean(rho_transfer_trials);
    sample_sz = size(rho_transfer_trials,1);
    ci_rho_transfer_trial = [std(rho_transfer_trials)./sqrt(sample_sz)] .* tinv(0.975,sample_sz-1); 
    shadedErrorBar([1:5],m_rho_transfer_trial,ci_rho_transfer_trial);
    xlim([0,7])
    xticks(1:5)
    xlabel('trials');
    ylabel('mean Fisher correlation coeff')


    figure;
    m_rho_block_rep = squeeze(mean(rho_block_rep,'omitnan'));
    ci_rho_block_rep = squeeze(std(rho_block_rep,'omitnan'))./sqrt(squeeze(sum(~isnan(rho_block_rep))));
    for b = 1:size(m_rho_block_rep,1)
        errorbar(m_rho_block_rep(b,:),ci_rho_block_rep(b,:),'CapSize',0,'Color',[0.8,0.8,0.8]); hold all;
        plot(m_rho_block_rep(b,:),'Color',cspace(b,:),'LineWidth',2);
        xlim([0,7])
        xticks(1:6)
        xlabel('stim repeat');
        ylabel('mean Fisher correlation coeff')
    end


    % median split based on transfer error behavior
    % behavior_residual_calc;
    % median_TE = median(transferErrorM);
    % med_split_idx = transferErrorM>median_TE;
    % figure;
    % m_rho_rep_l = mean(rho_stimOld_rep(med_split_idx==0,:));
    % m_rho_rep_h = mean(rho_stimOld_rep(med_split_idx==1,:));
    % 
    % ci_rho_rep_l = [std(rho_stimOld_rep(med_split_idx==0,:))./sqrt(sum(med_split_idx==0))] .* norminv(0.975);
    % ci_rho_rep_h = [std(rho_stimOld_rep(med_split_idx==1,:))./sqrt(sum(med_split_idx==1))] .* norminv(0.975);
    % 
    % errorbar(m_rho_rep_l,ci_rho_rep_l,'CapSize',0,'Color','r');
    % hold all;plot(m_rho_rep_l,'r','LineWidth',2);
    % 
    % errorbar(m_rho_rep_h,ci_rho_rep_h,'CapSize',0,'Color','b');
    % hold all;plot(m_rho_rep_h,'b','LineWidth',2);
    % 
    % xlim([0,7])
    % xticks(1:6)
    % xlabel('stim repeat');
    % ylabel('mean Fisher correlation coeff')



    % median split based on bt-state coefficients in another brain area
    btstate_coeff = ppc_btstate.gmean_b(:,1);
    median_btstate = median(btstate_coeff);
    med_btstate_split_idx = btstate_coeff>median_btstate;
    figure;
    m_rho_trial_l = mean(rho_stimOld_rep(med_btstate_split_idx==0,:));
    m_rho_trial_h = mean(rho_stimOld_rep(med_btstate_split_idx==1,:));

    ci_rho_rep_l = [std(rho_stimOld_rep(med_btstate_split_idx==0,:))./sqrt(sum(med_btstate_split_idx==0))] .* norminv(0.975);
    ci_rho_rep_h = [std(rho_stimOld_rep(med_btstate_split_idx==1,:))./sqrt(sum(med_btstate_split_idx==1))] .* norminv(0.975);
    
    errorbar(m_rho_trial_l,ci_rho_rep_l,'CapSize',0,'Color','r');
    hold all;plot(m_rho_trial_l,'r','LineWidth',2);

    errorbar(m_rho_trial_h,ci_rho_rep_h,'CapSize',0,'Color','b');
    hold all;plot(m_rho_trial_h,'b','LineWidth',2);

    xlim([0,7])
    xticks(1:6)
    xlabel('stim repeat');
    ylabel('mean Fisher correlation coeff')
    title('median split by btstate coeff')


%   split by bt state coefficients for trials in first miniblock

     figure;
    m_rho_trial_l = mean(rho_transfer_trials(med_btstate_split_idx==0,:));
    m_rho_trial_h = mean(rho_transfer_trials(med_btstate_split_idx==1,:));

    ci_rho_rep_l = [std(rho_transfer_trials(med_btstate_split_idx==0,:))./sqrt(sum(med_btstate_split_idx==0))] .* norminv(0.975);
    ci_rho_rep_h = [std(rho_transfer_trials(med_btstate_split_idx==1,:))./sqrt(sum(med_btstate_split_idx==1))] .* norminv(0.975);
    
    errorbar(m_rho_trial_l,ci_rho_rep_l,'CapSize',0,'Color','r');
    hold all;plot(m_rho_trial_l,'r','LineWidth',2);

    errorbar(m_rho_trial_h,ci_rho_rep_h,'CapSize',0,'Color','b');
    hold all;plot(m_rho_trial_h,'b','LineWidth',2);

    xlim([0,6])
    xticks(1:5)
    xlabel('trial');
    ylabel('mean Fisher correlation coeff')
    title('median split by btstate coeff')

    % split between correct in current and alternate state
    figure;
    rho_corrOld_rep_m = mean(rho_corrOld_rep,'omitnan');
    rho_corrAltOld_rep_m = mean(rho_corrAltOld_rep,'omitnan');
    rho_errorOld_rep_m = mean(rho_errorOld_rep,'omitnan');
    
    ci_rho_rep_c = [std(rho_corrOld_rep,'omitnan')./sqrt(sum(~isnan(rho_corrOld_rep)))] .* norminv(0.975);
    ci_rho_rep_ca = [std(rho_corrAltOld_rep,'omitnan')./sqrt(sum(~isnan(rho_corrAltOld_rep)))] .* norminv(0.975);
    ci_rho_rep_e = [std(rho_errorOld_rep,'omitnan')./sqrt(sum(~isnan(rho_errorOld_rep)))] .* norminv(0.975);

    errorbar(rho_corrOld_rep_m,ci_rho_rep_c,'CapSize',0,'Color','g');
    hold all;plot(rho_corrOld_rep_m,'g','LineWidth',2);

    errorbar(rho_corrAltOld_rep_m,ci_rho_rep_ca,'CapSize',0,'Color','b');
    hold all;plot(rho_corrAltOld_rep_m,'b','LineWidth',2);


    % 
    errorbar(rho_errorOld_rep_m,ci_rho_rep_e,'CapSize',0,'Color','r');
    hold all;plot(rho_errorOld_rep_m,'r','LineWidth',2);

    xlim([0,3])
    xticks(1:2)
    xlabel('stim repeat');
    ylabel('mean Fisher correlation coeff')

    title(sprintf('angle thresh = %d deg',angle_thresh))
    % legend({'corr','corrAlt','error'})

    figure;
    errorbar(mean(rho_window,'omitnan'),norminv(0.975)*(std(rho_window,'omitnan')./sqrt(size(rho_window,1))));
    xticks(1:length(window));
    xticklabels(window);
    xlim([0,length(window)+1]);
    hold all;plot([0,length(window)+1],[0,0],'r:');
    xlabel('error')
    ylabel('Fisher Correl Coeff');

    figure;
    errorbar(mean(rho_window_ref,'omitnan'),norminv(0.975)*(std(rho_window_ref,'omitnan')./sqrt(size(rho_window,1))));
    xticks(1:length(windowSigned));
    xticklabels(windowSigned);
    xlim([0,length(windowSigned)+1]);
    hold all;plot([0,length(windowSigned)+1],[0,0],'r:');
    xlabel('error')
    ylabel('Fisher Correl Coeff');


    rho0_window_signed_m = squeeze(mean(rho_window_signed(:,1,:),'omitnan'));
    rho0_window_signed_ci = norminv(0.975).*(squeeze(std(rho_window_signed(:,1,:),'omitnan'))./squeeze(sqrt(sum(~isnan(rho_window_signed(:,1,:))))));
    rho1_window_signed_m = squeeze(mean(rho_window_signed(:,2,:),'omitnan'));
    rho1_window_signed_ci = norminv(0.975)*(squeeze(std(rho_window_signed(:,2,:),'omitnan'))./squeeze(sqrt(sum(~isnan(rho_window_signed(:,2,:))))));

    figure;
    errorbar(rho0_window_signed_m,rho0_window_signed_ci,'b');hold all;
    errorbar(rho1_window_signed_m,rho1_window_signed_ci,'r'),
    xticks(1:length(windowSigned));
    xticklabels(windowSigned);
    xlim([0,length(windowSigned)+1]);
    plot([0,length(windowSigned)+1],[0,0],'r:');
    xlabel('error')
    ylabel('Fisher Correl Coeff');

end

disp('jee')

end

function fish = fisher_transform(data)
    fish = 0.5*log((1+data)./(1-data));
end

function betas = weighted_regress(Y,X,W)
    betas = inv(X'*W*X)*X'*W*Y;
end

