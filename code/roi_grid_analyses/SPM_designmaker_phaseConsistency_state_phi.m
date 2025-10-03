function SPM_designmaker_phaseConsistency_state_phi(subj, run, phi, roi, saveDir, numFolds)
% this function takes the events structure from gridplanet and
% converts it into a format that can be used as conditions for fMRI
% analysis with ~~ separate runs~

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


% sec/TR
TR=2;


[T, subData.store] = gridAnalysisFunction(subj,[1066,1080]);
T.TrialNum = [1:height(T)]';

XCenter = subData.xCenter;
YCenter = subData.yCenter;


%calculated from center of circle
idealAngle = atan2d((T.MeanY-YCenter),(T.MeanX-XCenter));
actualAngle = atan2d((T.RespY-YCenter),(T.RespX-XCenter));
idealAngle  = mod(idealAngle+360,360);
actualAngle  = mod(actualAngle+360,360);

T.angleDiff = abs(rad2deg(circ_dist(deg2rad(idealAngle),deg2rad(actualAngle))));
lenScanBlk = 84;

%      calculate distance travelled
T.travelDist = sqrt(((T.ShieldStartX-T.RespX).^2)+(T.ShieldStartY-T.ShieldStartY-T.RespY).^2);

T.errorExclude = T.angleDiff>45 ;%include only trials were degree of error is less than 45 (half of 90)
T.noMovement = T.travelDist < 10;

%for dropping trials analysis:
state2 = zeros(lenScanBlk,1);
state2(lenScanBlk/2+1:end) = 1;
if subj == 1066 || subj == 1080
    state=repmat(state2,3,1);
else
    state=repmat(state2,4,1);
end
T.state = state;




% separate out trials for current run
runT = T(T.MRIBlock == run,:);


% get events for current run
% run_events = subData.store(run).Event;
run_time = subData.store(run).Time;

% get indicies for sim, response and feedback in events list
sim_event_idx = find(contains(subData.store(run).Event,'sim_onset'));
resp_event_idx = find(contains(subData.store(run).Event,'resp_onset'));
outcome_event_idx = find(contains(subData.store(run).Event,'outcome_onset'));

% specify indicies for each trial type...
state_0_include = (runT.state == 0 & runT.errorExclude == 0 & runT.noMovement ==0);
state_0_exclude = (runT.state == 0 & runT.errorExclude == 1 & runT.noMovement ==0);
state_1_include = (runT.state == 1 & runT.errorExclude == 0 & runT.noMovement ==0);
state_1_exclude = (runT.state == 1 & runT.errorExclude == 1 & runT.noMovement ==0);
noMovement = runT.noMovement == 1;

% make parametric modulators for different conditions

modState0IncludeSpatial(:,1) = cos(numFolds*(deg2rad(runT.respAngle(state_0_include))-phi(2)));

modState1IncludeSpatial(:,1) = cos(numFolds*(deg2rad(runT.respAngle(state_1_include))-phi(1)));

modState0ExcludeSpatial(:,1) = cos(numFolds*(deg2rad(runT.respAngle(state_0_exclude))-phi(2)));

modState1ExcludeSpatial(:,1) = cos(numFolds*(deg2rad(runT.respAngle(state_1_exclude))-phi(1)));


% mean center
modState0IncludeSpatial = modState0IncludeSpatial-nanmean(modState0IncludeSpatial);
modState0IncludeSpatial(~isfinite(modState0IncludeSpatial))=0;

modState1IncludeSpatial = modState1IncludeSpatial-nanmean(modState1IncludeSpatial);
modState1IncludeSpatial(~isfinite(modState1IncludeSpatial))=0;

modState0ExcludeSpatial = modState0ExcludeSpatial-nanmean(modState0ExcludeSpatial);
modState0ExcludeSpatial(~isfinite(modState0ExcludeSpatial))=0;

modState1ExcludeSpatial = modState1ExcludeSpatial-nanmean(modState1ExcludeSpatial);
modState1ExcludeSpatial(~isfinite(modState1ExcludeSpatial))=0;



%    get times and convert to TRs
TR=2;
TimeTR = (run_time-run_time(1))/TR;


%   fill onsets and durations for conditions

%  ### SIMULATION PERIOD ##
%       state 0 - include
names{end+1} = 'state0_include';
onsets(end+1) =  {TimeTR(sim_event_idx(state_0_include))};
durations{end+1} = runT.WaitSimTime(state_0_include)./TR;

%       state 1 - include
names{end+1} = 'state1_include';
onsets(end+1) =  {TimeTR(sim_event_idx(state_1_include))};
durations{end+1} = runT.WaitSimTime(state_1_include)./TR;

%       state 0 - exclude
names{end+1} = 'state0_exclude';
onsets(end+1) =  {TimeTR(sim_event_idx(state_0_exclude))};
durations{end+1} = runT.WaitSimTime(state_0_exclude)./TR;

%       state 0 - exclude
names{end+1} = 'state1_exclude';
onsets(end+1) =  {TimeTR(sim_event_idx(state_1_exclude))};
durations{end+1} = runT.WaitSimTime(state_1_exclude)./TR;


%       no movement trials
names{end+1} = 'nomove';
onsets(end+1) =  {TimeTR(sim_event_idx(noMovement))};
durations{end+1} = runT.WaitSimTime(noMovement)./TR;

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
pmod(1).name{1}='State0IncSpa';
pmod(1).param{1}=modState0IncludeSpatial;
pmod(1).poly{1}=1;
%
pmod(2).name{1}='State1IncSpa';
pmod(2).param{1}=modState1IncludeSpatial;
pmod(2).poly{1}=1;
%
pmod(3).name{1}='State0ExcSpa';
pmod(3).param{1}=modState0ExcludeSpatial;
pmod(3).poly{1}=1;
%
pmod(4).name{1}='State1ExcSpa';
pmod(4).param{1}=modState1ExcludeSpatial;
pmod(4).poly{1}=1;
%

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
if any(empty_idx < 5)
    empty_idx_mod=empty_idx(any(empty_idx < 5));
    pmod(empty_idx_mod)=[];
end



%         set orthogonalization for regressors
orth=[repmat({false},1,length(names))];



%        set condition for orthogonalisation
filename=[saveDir sprintf('%d_gridplanet_phaseConsistency_state_run%d_%s.mat',subj,run,roi)];
save(filename,'names','onsets','durations','pmod','orth');