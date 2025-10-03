function SPM_designmaker_phaseConsistency_state_allruns_phi(subj, phi, roi, saveDir, numFolds)
% this function takes the events structure from gridplanet and
% converts it into a format that can be used as conditions for fMRI
% analysis with ~~ all runs concatenated ~~

% files directory
filedir = sprintf('/gpfs/data/mnassar/lyu21/GridPlanet/');

addpath(filedir)

% load raw behavioral data
if subj == 1066 || subj == 1080
    subData=load([filedir sprintf('behavioralData/gridplanet_Scanning_subj_%d_block_3.mat',subj)]);
else
    subData=load([filedir sprintf('behavioralData/gridplanet_Scanning_subj_%d_block_4.mat',subj)]);

end


% make variables for conditions
durations={};
onsets={};
names={};
% make empty structure for parametric modulator
pmod = struct('name',(''),'param',{},'poly',{});


% seconds/TR
TR=2;


[T, subData.store] = gridAnalysisFunction(subj,[1066,1080]);
T.TrialNum = [1:height(T)]';

XCenter = subData.xCenter;
YCenter = subData.yCenter;
  
    
    
    %get error
%      distMean =  sqrt((T.MeanX-T.RespX).^2 + (T.MeanY - T.RespY).^2);    
    %calculated from center of circle
    idealAngle = atan2d((T.MeanY-YCenter),(T.MeanX-XCenter));
     actualAngle = atan2d((T.RespY-YCenter),(T.RespX-XCenter));
      idealAngle  = mod(idealAngle+360,360);
      actualAngle  = mod(actualAngle+360,360);
     
   T.angleDiff = abs(rad2deg(circ_dist(deg2rad(idealAngle),deg2rad(actualAngle))));
     lenScanBlk = 84;
     T.idealAngle = idealAngle;
     
     %     calculate error for alternative state
    idealAngleAlt = zeros(length(idealAngle),1);
    for k = 1:height(T)
%         %     calculate angle to alternative state for each trial
%         idealAngleAlt(k,1) = atan2d(color_state_meanY(T.ColorID(k),states~=T.state(k))-YCenter,...
%             color_state_meanX(T.ColorID(k),states~=T.state(k))-XCenter);
%         idealAngleAlt(k,1)  = mod(idealAngleAlt(k,1)+360,360);
         if T.Block(k) == 1
             idealAngleAlt(k,1) = mod(T.idealAngle(k) - 90,360);
         else
             idealAngleAlt(k,1) = mod(T.idealAngle(k) + 90,360);
         end
    end

    T.angleDiffAlt = abs(rad2deg(circ_dist(deg2rad(idealAngleAlt),deg2rad(actualAngle))));
     
%      calculate distance travelled 
    T.travelDist = sqrt(((T.ShieldStartX-T.RespX).^2)+(T.ShieldStartY-T.RespY).^2);


%    trials where angular distance is not in line with either state
    angleThresh = 45;
    moveThresh = 10;
    T.noMovement = T.travelDist <= moveThresh;
    T.errorExclude = T.angleDiff>angleThresh & T.angleDiffAlt > angleThresh & T.noMovement == 0;  
    T.errorExclude = T.angleDiff>angleThresh & T.angleDiffAlt > angleThresh & T.noMovement == 0;
    T.correctAngle = T.angleDiff<angleThresh & T.angleDiffAlt >=angleThresh & T.noMovement == 0;
    T.correctAngleAlt = T.angleDiff>=angleThresh & T.angleDiffAlt<angleThresh & T.noMovement == 0;

     %for dropping trials analysis:
            state2 = zeros(lenScanBlk,1);
            state2(lenScanBlk/2+1:end) = 1;
            if subj == 1066 || subj == 1080
                state=repmat(state2,3,1);
            else
                state=repmat(state2,4,1);
            end
            T.state = state;
            

if subj == 1066 || subj == 1080
    numruns = 1:3;
else
    numruns = 1:4;
end



% concatenate run events and their timing
run_events = [];
run_time =[];
for run = numruns
    run_events = [run_events;subData.store(run).Event];
    run_time_temp = subData.store(run).Time;
%     subtract block start time from each run
    run_time_temp = run_time_temp-run_time_temp(1);
%     add in the number of TRs for previous runs
    run_time_temp = run_time_temp + 346*(run-1)*TR;
    run_time = [run_time;run_time_temp];
end


% get indicies for sim, response and feedback in events list
sim_event_idx = find(contains(run_events,'sim_onset'));
resp_event_idx = find(contains(run_events,'resp_onset'));
outcome_event_idx = find(contains(run_events,'outcome_onset'));


% specify indicies for each trial type...
state_0_corr = (T.state == 0 & T.correctAngle == 1 & T.noMovement == 0);
state_1_corr = (T.state == 1 & T.correctAngle == 1 & T.noMovement == 0);
state_0_error = (T.state==0 & T.errorExclude == 1 & T.noMovement == 0);
state_1_error = (T.state==1 & T.errorExclude == 1 & T.noMovement == 0);
noMovement = T.noMovement == 1;

% make parametric modulators for different conditions
modState0corrSpatial(:,1) = cos(numFolds*(deg2rad(T.respAngle(state_0_corr))-phi(2)));

modState1corrSpatial(:,1) = cos(numFolds*(deg2rad(T.respAngle(state_1_corr))-phi(1)));

modState0ErrorSpatial(:,1) = cos(numFolds*(deg2rad(T.respAngle(state_0_error))-phi(2)));

modState1ErrorSpatial(:,1) = cos(numFolds*(deg2rad(T.respAngle(state_1_error))-phi(1)));



% mean center
modState0corrSpatial = modState0corrSpatial-nanmean(modState0corrSpatial);
modState0corrSpatial(~isfinite(modState0corrSpatial))=0;

modState1corrSpatial = modState1corrSpatial-nanmean(modState1corrSpatial);
modState1corrSpatial(~isfinite(modState1corrSpatial))=0;

modState0ErrorSpatial = modState0ErrorSpatial-nanmean(modState0ErrorSpatial);
modState0ErrorSpatial(~isfinite(modState0ErrorSpatial))=0;

modState1ErrorSpatial = modState1ErrorSpatial-nanmean(modState1ErrorSpatial);
modState1ErrorSpatial(~isfinite(modState1ErrorSpatial))=0;


%    get times and convert to TRs
% TR=2;
TimeTR = (run_time-run_time(1))/TR;

    
%   fill onsets and durations for conditions

%  ### SIMULATION PERIOD ##
%       state 0 - include
        names{end+1} = 'state0_corr';
        onsets(end+1) =  {TimeTR(sim_event_idx(state_0_corr))};
        durations{end+1} = T.WaitSimTime(state_0_corr)./TR;

%       state 1 - include
        names{end+1} = 'state1_corr';
        onsets(end+1) =  {TimeTR(sim_event_idx(state_1_corr))};
        durations{end+1} = T.WaitSimTime(state_1_corr)./TR;


% %      errors
        if sum(state_0_error) >  0
            names{end+1} = 'state_0_error';
            onsets(end+1) =  {TimeTR(sim_event_idx(state_0_error))};
            durations{end+1} = T.WaitSimTime(state_0_error)./TR;
        end
        if sum(state_1_error) >  0
            names{end+1} = 'state_1_error';
            onsets(end+1) =  {TimeTR(sim_event_idx(state_1_error))};
            durations{end+1} = T.WaitSimTime(state_1_error)./TR;
        end
        if sum(noMovement) >  0
            names{end+1} = 'nomove';
            onsets(end+1) =  {TimeTR(sim_event_idx(noMovement))};
            durations{end+1} = T.WaitSimTime(noMovement)./TR;
        end
        
%         #### RESPONSES

%         responses 
        names{end+1} = 'resp';
        onsets(end+1) =  {TimeTR(resp_event_idx)};
        durations{end+1} = repmat(3/TR,length(onsets{end}),1);        

%          #### FEEDBACK
%         feedback 
        names{end+1} = 'fdbk';
        onsets(end+1) =  {TimeTR(outcome_event_idx)};
        durations{end+1} = repmat(1/TR,length(onsets{end}),1);    

        
        
%         add parametric modulators
       pmod(1).name{1}='State0CorrSpa';
       pmod(1).param{1}=modState0corrSpatial;
       pmod(1).poly{1}=1;
% 
       pmod(2).name{1}='State1CorrSpa';
       pmod(2).param{1}=modState1corrSpatial;
       pmod(2).poly{1}=1;
% 
       if sum(state_0_error) > 0
           pmod(3).name{1}='State0ErrorSpa';
           pmod(3).param{1}=modState0ErrorSpatial;
           pmod(3).poly{1}=1;
       end
% 
       if sum(state_1_error) > 0
           pmod(4).name{1}='State1ErrorSpa';
           pmod(4).param{1}=modState1ErrorSpatial;
           pmod(4).poly{1}=1;
       end 

%         loop through regressors and remove any without any trials....
        empty_idx = [];
        for j = 1:length(onsets)
            if isempty(onsets{j})
                empty_idx(end+1) = j;
            end
        end
%             remove empty regressors
        onsets(empty_idx)=[];
        durations(empty_idx)=[];
        names(empty_idx)=[];
        pmod(empty_idx)=[];
        
        
        
%         set orthogonalization for regressors
        orth=[repmat({false},1,length(names))];



%        set condition for orthogonalisation
if numFolds == 6
    filename=[saveDir sprintf('%d_gridplanet_phaseConsistency_state_allruns_%s.mat',subj,roi)];
else
    filename=[saveDir sprintf('%d_gridplanet_phaseConsistency_state_allruns_%s_%dFold.mat',subj,roi,numFolds)];
end
save(filename,'names','onsets','durations','pmod','orth'); 