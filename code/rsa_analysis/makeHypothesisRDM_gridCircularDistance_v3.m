
function [R, exclude] = makeHypothesisRDM_gridCircularDistance_v3(T,subData, trialSelect, excludeError, excludeRep)
%this function makes hypothesis regressors for grid planet fmri study


numColors=7;
numStimReps=max(unique(T.StimulusRep));
numtrials = height(T);
%

if subData.subjNum == 1066 || subData.subjNum == 1080
    numruns = 3;
else
    numruns = 4;
end

trialsperState = numColors*numStimReps;

if numruns == 4
    T.stateVec = [ones(trialsperState,1);ones(trialsperState,1)*2; ...
        ones(trialsperState,1);ones(trialsperState,1)*2;...
        ones(trialsperState,1);ones(trialsperState,1)*2; ...
        ones(trialsperState,1);ones(trialsperState,1)*2];
elseif numruns == 3
    T.stateVec = [ones(trialsperState,1);ones(trialsperState,1)*2; ...
        ones(trialsperState,1);ones(trialsperState,1)*2;...
        ones(trialsperState,1);ones(trialsperState,1)*2];
end

%4. scanner run (need to create) 80 trials per run, 6 runs
lengthRun=sum(T.MRIBlock==1);
if numruns == 4
    T.scannerRuns = [ones(lengthRun,1);ones(lengthRun,1)*2;ones(lengthRun,1)*3;ones(lengthRun,1)*4];
elseif numruns == 3
    T.scannerRuns = [ones(lengthRun,1);ones(lengthRun,1)*2;ones(lengthRun,1)*3];
end

%amount of movement
%this takes the direction vector (at which timepoints is the participant
%moving vs. not moving) and averages them across the 5 trials in a miniblock
%motor movement per trial
direction = subData.direction;
T.amountofMovement =  sum(direction~=0)';


XCenter = subData.xCenter;
YCenter = subData.yCenter;



%get error
distMean =  sqrt((T.MeanX-T.RespX).^2 + (T.MeanY - T.RespY).^2);
%calculated from center of circle
idealAngle = atan2d((T.MeanY-YCenter),(T.MeanX-XCenter));
actualAngle = atan2d((T.RespY-YCenter),(T.RespX-XCenter));
idealAngle  = mod(idealAngle+360,360);
actualAngle  = mod(actualAngle+360,360);



% record state for first and second presentation of each novel color so
% that we can know its initial state and transfer state
NovelColorStateOrder = nan(13,2);
% make a vector for noting transition trials
transition_trial=zeros(height(T),1);
for j = 1:height(T)
    %     identify novel colors and save their initial and transfer states...
    if T.ColorID(j) > 5
        if isnan(NovelColorStateOrder(T.ColorID(j),1))
            NovelColorStateOrder(T.ColorID(j),1)=T.stateVec(j);
        elseif T.stateVec(j) ~=  NovelColorStateOrder(T.ColorID(j),1)
            NovelColorStateOrder(T.ColorID(j),2) = T.stateVec(j);
        end
    end
    % identify transfer trials
    if j > 1
        if T.stateVec(j) ~= T.stateVec(j-1)
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
    state_change_trial(transition_idx(j):transition_idx(j)+6)=1;
end

%     determine whether novel color trial is new or transfer trial
for j = 1:height(T)
    if T.ColorID(j) > 5
        if T.stateVec(j) == NovelColorStateOrder(T.ColorID(j),2) && state_change_trial(j) == 1
            transfer_trial(j) = 1;
        end
    end
end


exclude = zeros(height(T),1);
T.exclude = zeros(height(T),1);


%get errors
distMean =  sqrt((T.MeanX-T.RespX).^2 + (T.MeanY - T.RespY).^2);

%calculated from center of circle
idealAngle = atan2d((T.MeanY-YCenter),(T.MeanX-XCenter));
actualAngle = atan2d((T.RespY-YCenter),(T.RespX-XCenter));
idealAngle  = mod(idealAngle+360,360);
actualAngle  = mod(actualAngle+360,360);
T.angleDiff = abs(rad2deg(circ_dist(deg2rad(idealAngle),deg2rad(actualAngle))));


if excludeError ==1
    T.exclude = T.angleDiff>45; %include only trials were degree of error is less than 45 (half of 90)
    %remove trials where there were no movement
    for i = 1:height(T)
        if T.ShieldStartX(i) == T.RespX(i) && T.ShieldStartY(i) == T. RespY(i)
            %exclude(i)=1;
            T.exclude(i) = 1;
        end
    end
    
    exclude = T.exclude;
    T((T.exclude==1),:)=[];
    
end

if excludeRep > 0
    T.exclude(T.StimulusRep > excludeRep) = 1;
    exclude = T.exclude;
    T(T.exclude==1,:)=[];
end



%first, get the trajectory in radians.
trajAngle = round(deg2rad(T.respAngle),2);
%cognitive (adding state 2 vecs to 90)
cogAngle = T.respAngle;
cogAngle(T.stateVec==2) = cogAngle(T.stateVec==2)+90;
cogAngle(cogAngle>360)=cogAngle(cogAngle>360)-360;
cogAngle = round(deg2rad(cogAngle),2);

%create index.
a = round([0:0.01:2*pi],2);

repLim = 3;

%color by state
for i = 1:height(T)
    for j = 1:height(T)
        
        if T.ColorID(i) == T.ColorID(j) && mod(T.Block(i),2)~=mod(T.Block(j),2) && T.StimulusRep(i) > repLim && T.StimulusRep(j) <= repLim
            colorBetweenStateEarlyRDM(i,j) = -1;
        else
            colorBetweenStateEarlyRDM(i,j) = 0;
        end
        
        if T.ColorID(i) == T.ColorID(j) && mod(T.Block(i),2)==mod(T.Block(j),2) %same color, same state
            colorStateRDM(i,j)=-1;
            colorWithinStateRDM(i,j) = -1;
            colorBetweenStateRDM(i,j) = 0;
        elseif T.ColorID(i) == T.ColorID(j) && mod(T.Block(i),2)~=mod(T.Block(j),2) %if they are in the same color but in different states
            colorStateRDM(i,j)=-1;
            colorWithinStateRDM(i,j) = 0;
            colorBetweenStateRDM(i,j) = -1;
        elseif T.ColorID(i) ~= T.ColorID(j) && mod(T.Block(i),2)==mod(T.Block(j),2)  %if they are different colors, in the same state
            colorStateRDM(i,j)=1;
            colorWithinStateRDM(i,j) = 1;
            colorBetweenStateRDM(i,j) = 0;
        elseif T.ColorID(i) ~= T.ColorID(j) && mod(T.Block(i),2)~=mod(T.Block(j),2)  %if they are different colors, in different states
            colorStateRDM(i,j)=1;
            colorWithinStateRDM(i,j) = 0;
            colorBetweenStateRDM(i,j) = 1;
        end
    end
end



%color only
for i = 1:height(T)
    for j = 1:height(T)
        if T.ColorID(i) == T.ColorID(j)
            colorRDM(i,j)=0;
        else
            colorRDM(i,j)=1;
        end
    end
end

allTrials = [T.Block, T.stateVec, T.scannerRuns];
block=1;stateS=2;scanR=3;

numTrials = height(T);

%initialize RDM
%cueRDM = zeros(numTrials,numTrials);
runRDM = zeros(numTrials,numTrials);
blockRDM = zeros(numTrials,numTrials);
stateRDM = zeros(numTrials,numTrials);
movementRDM=zeros(numTrials,numTrials);


%make motor movement RSA
%now make a dissimilarity (distance matrix), by taking pairwise distances
%with itself (could also do m-m')
for m = 1:numTrials
    for n=1:numTrials
        movementRDM(m,n) = abs(T.amountofMovement(m)-T.amountofMovement(n));
    end
end

T.MovementX=abs(T.RespX-T.ShieldStartX);
T.MovementY=abs(T.RespY-T.ShieldStartY);

% get 'invariant' positions for the shields...
colorIDs = unique(T.ColorID);
cx0 = zeros(length(colorIDs),1);
cy0 = zeros(length(colorIDs),1);
ca0 = zeros(length(colorIDs),1);
cx1 = zeros(length(colorIDs),1);
cy1 = zeros(length(colorIDs),1);
ca1 = zeros(length(colorIDs),1);


% convert shield colors to CIE Lab space
colorRGB = [subData.shieldColors;subData.addShieldColors];
colorLAB = rgb2lab(colorRGB);

for c = 1:length(colorIDs)
    tempX0 = T.MeanX(T.ColorID==colorIDs(c) & mod(T.Block,2) ~= 0);
    tempY0 = T.MeanY(T.ColorID==colorIDs(c) & mod(T.Block,2) ~= 0);
    tempA0 = idealAngle(T.ColorID==colorIDs(c) & mod(T.Block,2) ~= 0);
    tempX1 = T.MeanX(T.ColorID==colorIDs(c) & mod(T.Block,2) == 0);
    tempY1 = T.MeanY(T.ColorID==colorIDs(c) & mod(T.Block,2) == 0);
    tempA1 = idealAngle(T.ColorID==colorIDs(c) & mod(T.Block,2) == 0);
    
    % if we cannot find a presentation of this color in first state, then
    % use position from second state rotated back into first state
    % position... for subjs 1066, 1080
    if isempty(tempX0)
        tempX1 = T.MeanX(T.ColorID==colorIDs(c) & mod(T.Block,2) == 0)-XCenter;
        tempY1 = T.MeanY(T.ColorID==colorIDs(c) & mod(T.Block,2) == 0)-YCenter;
        tempA1 = idealAngle(T.ColorID==colorIDs(c) & mod(T.Block,2) == 0);
        tempX0 = (tempX1(1)*cos(pi/2) + tempY1(1)*sin(pi/2))+XCenter;
        tempY0 = -(tempX1(1)*sin(pi/2) + tempY1(1)*cos(pi/2))+YCenter;
        tempA0 = mod(tempA1(1) - 90,360);
    end
    if isempty(tempX1)
        tempX0 = T.MeanX(T.ColorID==colorIDs(c) & mod(T.Block,2) ~= 0)-XCenter;
        tempY0 = T.MeanY(T.ColorID==colorIDs(c) & mod(T.Block,2) ~= 0)-YCenter;
        tempA0 = idealAngle(T.ColorID==colorIDs(c) & mod(T.Block,2) ~= 0);
        tempX1 = (tempX0(1)*cos(pi/2) + tempY0(1)*sin(pi/2))+XCenter;
        tempY1 = -(tempX0(1)*sin(pi/2) + tempY0(1)*cos(pi/2))+YCenter;
        tempA1 = mod(tempA0(1) - 90,360);
    end
    
    cx0(c,1) = tempX0(1);
    cy0(c,1) = tempY0(1);
    ca0(c,1) = tempA0(1);
    cx1(c,1) = tempX1(1);
    cy1(c,1) = tempY1(1);
    ca1(c,1) = tempA1(1);
end

% make a vector for alternate state positions for each trial...
for t = 1:numTrials
    if mod(T.Block(t),2) ~= 0
        meanErrorX(t,1) = cx1((T.ColorID(t)==colorIDs));
        meanErrorY(t,1) = cy1((T.ColorID(t)==colorIDs));
        errorAngle(t,1) = ca1(T.ColorID(t)==colorIDs);
    elseif mod(T.Block(t),2) == 0
        meanErrorX(t,1) = cx0((T.ColorID(t)==colorIDs));
        meanErrorY(t,1) = cy0((T.ColorID(t)==colorIDs));
        errorAngle(t,1) = ca0((T.ColorID(t)==colorIDs));
    end
end

%response and outcome
for m = 1:numTrials
    for n=1:numTrials
        startShield_RDM(m,n) =  sqrt((T.ShieldStartX(m)-T.ShieldStartX(n)).^2 + (T.ShieldStartY(m) - T.ShieldStartY(n)).^2);
        endShield_RDM(m,n) =  sqrt((T.RespX(m)-T.RespX(n)).^2 + (T.RespY(m) - T.RespY(n)).^2);
        movementShield_RDM(m,n) =  sqrt((T.MovementX(m)-T.MovementX(n)).^2 + (T.MovementY(m) - T.MovementY(n)).^2);
        response_RDM(m,n) =  sqrt((T.RespX(m)-T.RespX(n)).^2 + (T.RespY(m) - T.RespY(n)).^2);
        outcome_RDM(m,n) = sqrt((T.OutcomeX(m)-T.OutcomeX(n)).^2 + (T.OutcomeY(m) - T.OutcomeY(n)).^2);
        spatialMeanPosition_RDM(m,n) = sqrt((T.MeanX(m)-T.MeanX(n)).^2 + (T.MeanY(m) - T.MeanY(n)).^2);
        
        cogMeanPosition_RDM(m,n) = sqrt((cx0((T.ColorID(m)==colorIDs),1)-cx0((T.ColorID(n)==colorIDs),1)).^2 + ...
            (cy0((T.ColorID(m)==colorIDs),1)-cy0((T.ColorID(n)==colorIDs),1)).^2);
        
        colorLAB_RDM(m,n) = sqrt((colorLAB(T.ColorID(m)==colorIDs,1)-colorLAB(T.ColorID(n)==colorIDs,1)).^2 + ...
            (colorLAB(T.ColorID(m)==colorIDs,2)-colorLAB(T.ColorID(n)==colorIDs,2)).^2 + ...
            (colorLAB(T.ColorID(m)==colorIDs,3)-colorLAB(T.ColorID(n)==colorIDs,3)).^2);
        
        errorPosition_RDM(m,n) = sqrt((meanErrorX(m) - meanErrorX(n)).^2 + (meanErrorY(m) - meanErrorY(n)).^2);
        
        errorAngle_RDM(m,n) = mod(circ_dist(deg2rad(errorAngle(m)), deg2rad(errorAngle(n))),2*pi);
        spatialAngle_RDM(m,n) = mod(circ_dist(deg2rad(idealAngle(m)),deg2rad(idealAngle(n))),2*pi);
        cogAngle_RDM(m,n) = mod(circ_dist(deg2rad(ca0((T.ColorID(m)==colorIDs),1)),deg2rad(ca0((T.ColorID(n)==colorIDs),1))),2*pi);
    end
end




%DISSIMILARITY RDM = higher is less similar
%if using euclidean distance would want to do dissimilarity)
%if not a match == 2; if a match == 1;
for m = 1:length(allTrials)
    for n=1:length(allTrials)
        if allTrials(m, scanR)==allTrials(n,scanR)
            runRDM(m,n)= 1;
        else
            runRDM(m,n)=2;
        end
        if allTrials(m, block)==allTrials(n,block)
            blockRDM(m,n)= 1;
        else
            blockRDM(m,n)=2;
        end
        if allTrials(m, stateS)==allTrials(n,stateS)
            stateRDM(m,n)= 1;
        else
            stateRDM(m,n)=2;
        end
        
    end
end

% make regressors for transfer trials
transfer_trial_RDM = abs(transfer_trial-transfer_trial');


% make regressors for angle/pos changes
colorAngle=nan(13,1);
colorMeanX = nan(13,1);
colorMeanY = nan(13,1);
lastAngle = nan(13,1);

for i = 1:height(T)
    idealAngleChange(i) = circ_dist(deg2rad(idealAngle(i)),deg2rad(colorAngle(T.ColorID(i))));
    meanChange(i,1) = sqrt(((T.MeanX(i)-colorMeanX(T.ColorID(i)))^2) + ((T.MeanY(i)-colorMeanY(T.ColorID(i)))^2));
    angleChange(i) = circ_dist(deg2rad(actualAngle(i)),deg2rad(lastAngle(T.ColorID(i))));
    lastAngle(T.ColorID(i)) = actualAngle(i);
    colorAngle(T.ColorID(i)) = idealAngle(i);
    colorMeanX(T.ColorID(i)) = T.MeanX(i);
    colorMeanY(T.ColorID(i)) = T.MeanY(i);
end

angleChange(isnan(angleChange)) = 0;
meanChange(isnan(meanChange)) = 0;
meanChange_RDM = (abs(meanChange-meanChange'));
meanChangeSq_RDM = (abs(meanChange.^2-meanChange.^2'));


for i = 1:length(angleChange)
    for j = 1:length(angleChange)
        angleChange_RDM(i,j) = abs(circ_dist(angleChange(i),angleChange(j)));
    end
end


idealAngleChange(isnan(idealAngleChange))=0;
idealAngleChange_RDM = abs(idealAngleChange-idealAngleChange');
idealAngleSqChange_RDM = abs((idealAngleChange.^2) - (idealAngleChange.^2)');


%lagRDM=struct;

%create a lag RDM, that is n = n+1;
for l = 1:16 %number of lags I want
    lagRDM(:,:) = zeros(numtrials,numtrials); %initialize at full size
    for m = 1:length(lagRDM)
        for n=1:length(lagRDM)
            if m==n+l
                lagRDM(m,n)= 1;
            else
                lagRDM(m,n)=2;
            end
        end
    end
    %reduce size of RDM
    if sum(trialSelect) > 0
        lagum = lagRDM(trialSelect==1,:);
        lagRegressors(:,:,l) = reshape(lagum',[size(lagum,1)*size(lagum,2),1]);
    else
        lagum = lagRDM(exclude==0,exclude==0);
        %take lower triangle and make into regressor
        lagRegressors(:,:,l) = lagum(find(tril(lagum,-1)));
        
        
    end
end



%
%now take lower triangle (not include diagonal) make them into vectors
if sum(trialSelect) > 0
    runRDM = runRDM(trialSelect==1,:);
    blockRDM = blockRDM(trialSelect==1,:);
    stateRDM = stateRDM(trialSelect==1,:);
    movementRDM = movementRDM(trialSelect==1,:);
    movementShield_RDM = movementShield_RDM(trialSelect==1,:);
    colorStateRDM = colorStateRDM(trialSelect==1,:);
    colorBetweenStateRDM = colorBetweenStateRDM(trialSelect==1,:);
    colorRDM = colorRDM(trialSelect==1,:);
    startShield_RDM = startShield_RDM(trialSelect==1,:);
    outcome_RDM = outcome_RDM(trialSelect==1,:);
    endShield_RDM = endShield_RDM(trialSelect==1,:);
    
    R=struct;
    R.runRegressor = reshape(runRDM',[size(runRDM,1)*size(runRDM,2),1]);
    R.blockRegressor = reshape(blockRDM',[size(runRDM,1)*size(runRDM,2),1]);
    R.stateRegressor= reshape(stateRDM',[size(runRDM,1)*size(runRDM,2),1]);
    
    R.movementRegressor = reshape(movementRDM',[size(runRDM,1)*size(runRDM,2),1]);
    R.movementShieldRegressor = reshape(movementShield_RDM',[size(runRDM,1)*size(runRDM,2),1]);
    
    R.colorStateRegressor = reshape(colorStateRDM',[size(runRDM,1)*size(runRDM,2),1]);
    R.colorBetweenStateRegressor = reshape(colorBetweenStateRDM',[size(runRDM,1)*size(runRDM,2),1]);
    R.colorWithinStateRegressor = reshape(colorWithinStateRDM',[size(runRDM,1)*size(runRDM,2),1]);
    R.colorRegressor = reshape(colorRDM',[size(runRDM,1)*size(runRDM,2),1]);
    R.colorLAB_RDM = reshape(colorLAB_RDM',[size(runRDM,1)*size(runRDM,2),1]);
    
    R.startShieldRegressor = reshape(startShield_RDM',[size(runRDM,1)*size(runRDM,2),1]);
    R.outcomeRegressor = reshape(outcome_RDM',[size(runRDM,1)*size(runRDM,2),1]);
    R.endShieldRegressor = reshape(endShield_RDM',[size(runRDM,1)*size(runRDM,2),1]);
    R.spatialMeanPositionRegressor = reshape(spatialMeanPosition_RDM',[size(runRDM,1)*size(runRDM,2),1]);
    R.cogMeanPositionRegressor = reshape(cogMeanPosition_RDM',[size(runRDM,1)*size(runRDM,2),1]);
    R.errorPositionRegressor = reshape(errorPosition_RDM',[size(runRDM,1)*size(runRDM,2),1]);
    
    R.spatialAngleRegressor = reshape(spatialAngle_RDM',[size(runRDM,1)*size(runRDM,2),1]);
    R.cogAngleRegressor = reshape(cogAngle_RDM',[size(runRDM,1)*size(runRDM,2),1]);
    R.cogAngleRegressor = reshape(cogAngle_RDM',[size(runRDM,1)*size(runRDM,2),1]);
    R.errorAngleRegressor = reshape(errorAngle_RDM',[size(runRDM,1)*size(runRDM,2),1]);
    
    R.idealAngleChangeRegressor = reshape(idealAngleChange_RDM',[size(runRDM,1)*size(runRDM,2),1]);
    R.idealAngleSqChangeRegressor = reshape(idealAngleSqChange_RDM',[size(runRDM,1)*size(runRDM,2),1]);
    
    R.angleChangeRegressor = reshape(angleChange_RDM',[size(runRDM,1)*size(runRDM,2),1]);
    R.transferTrialRegressor = reshape(transfer_trial_RDM',[size(runRDM,1)*size(runRDM,2),1]);
    
    
else
    
    R=struct;
    R.runRegressor= runRDM(find(tril(runRDM,-1)));
    R.blockRegressor = blockRDM(find(tril(blockRDM,-1)));
    R.stateRegressor=stateRDM(find(tril(stateRDM,-1)));
    
    mask=tril(true(size(movementRDM)),-1);
    
    R.movementRegressor = movementRDM(mask);
    R.movementShieldRegressor = movementShield_RDM(mask);
    
    R.colorStateRegressor = colorStateRDM(mask);
    R.colorBetweenStateRegressor = colorBetweenStateRDM(mask);
    R.colorBetweenStateEarlyRegressor = colorBetweenStateEarlyRDM(mask);
    R.colorWithinStateRegressor = colorWithinStateRDM(mask);
    R.colorRegressor = colorRDM(mask);
    R.colorLABRegressor = colorLAB_RDM(mask);
    
    
    R.startShieldRegressor = startShield_RDM(mask);
    % R.responseRegressor=response_RDM(mask);
    R.outcomeRegressor = outcome_RDM(mask);
    R.endShieldRegressor = endShield_RDM(mask);
    R.spatialMeanPositionRegressor = spatialMeanPosition_RDM(mask);
    R.cogMeanPositionRegressor = cogMeanPosition_RDM(mask);
    R.errorPositionRegressor = errorPosition_RDM(mask);
    R.spatialAngleRegressor = spatialAngle_RDM(mask);
    R.cogAngleRegressor = cogAngle_RDM(mask);
    R.errorAngleRegressor = errorAngle_RDM(mask);
    R.idealAngleChangeRegressor = idealAngleChange_RDM(mask);
    R.idealAngleSqChangeRegressor = idealAngleSqChange_RDM(mask);
    R.angleChangeRegressor = angleChange_RDM(mask);
    R.transferTrialRegressor = transfer_trial_RDM(mask);
    
    
    
end


R.lagRegressors = lagRegressors;


end






