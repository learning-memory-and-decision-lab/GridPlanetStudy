% this script runs analyses on behavioral data to compare performance on
% transfer trials to random guessing

clear all

% ------- parameters for analysis - num iterations, session number, etc
% Analysis parameters
nIt = 1000;

% session number
session = 2;

% novel color setting (session 2)
novelColors = true;
% -------------------------------------------------------------

% Set path
if session == 1 || session == 3
    dataDir = sprintf('~/Dropbox/Gridplanet/Data/session%d/',session);
else
    dataDir = sprintf('~/Dropbox/Gridplanet/Data/session%d/extracted_data/',session);
end
cd(dataDir)

if session == 1 || session == 3
    a=dir('sortedData*')
else
    a=dir('*.mat')
end

% get file names:
fns={a.name};


% Loop through files and extract errors and trial labels:
for i =1:length(fns)
    fn=fns{i}
    subData=load(fn);
    if session == 1
        ss = fn(19:22);
        load(fullfile(dataDir,sprintf('rawData/gridplanet_V2_subj_%s.mat',ss)),'xCenter','yCenter');
    elseif session == 2
        ss = fn(1:4);
        load(fullfile(dataDir,sprintf('../rawData/gridplanet_Scanning_subj_%s_block_3.mat',ss)),'xCenter','yCenter');
    elseif session == 3
        ss = fn(19:22);
        load(fullfile(dataDir,sprintf('rawData/gridplanet_V2_subj_%s_angleswitch.mat',ss)),'xCenter','yCenter');
    end
    
    % Calculate angles:
    idealAngle = deg2rad(atan2d((subData.T.MeanY-yCenter),(subData.T.MeanX-xCenter)));
    actualAngle = deg2rad(atan2d((subData.T.RespY-yCenter),(subData.T.RespX-xCenter)));
    subData.T.idealAngle = idealAngle;
    subData.T.SN=repmat(i, length(subData.T.Block), 1);
    
    
%   if this is session 1 or 2, the angular shift between states is always
%   90 degrees. However, if this is session 3, then the angular shift is
%   changing on each block and we need to calculate 
    if session == 1 || session == 2
        angleShift = pi/2;
    elseif session == 3
        blockNums = unique(subData.T.Block);
        colorNums = unique(subData.T.ColorID);
        
        for b = blockNums'
            for c = colorNums'
                colorMeanX_by_Block(c,b) = mean(subData.T.MeanX(subData.T.ColorID==c & subData.T.Block==b));
                colorMeanY_by_Block(c,b) = mean(subData.T.MeanY(subData.T.ColorID==c & subData.T.Block==b));
            end
        end
%         now run through the series of blocks and calculate angular
%         distance from center
        for b = 1:length(blockNums)
            angleBlock(b,1) = deg2rad(mod(atan2d((colorMeanY_by_Block(1,b)-yCenter),(colorMeanX_by_Block(1,b)-xCenter)),360));
        end
%         now calculate the angular distance between blocks
        for b = 2:length(blockNums)
            angleChangeBlock(b,1) = (mod((circ_dist((angleBlock(b)),(angleBlock(b-1)))),2*pi));
        end
    end
                
            
        
    %     calculate error for alternative state
        idealAngleAlt = zeros(length(idealAngle),1);
        for k = 1:height(subData.T)
%             if this is session 3, then flexibly update angleShift by
%             block number
            if session == 3
                if subData.T.Block(k) == 1
                    angleShift = pi/2;
                else
                    angleShift = angleChangeBlock(subData.T.Block(k));
                end
            end

            if session == 1 || session == 2
                 if subData.T.Block(k) == 1 || subData.T.Block(k) == 3 || subData.T.Block(k) == 5 || subData.T.Block(k) == 7 
                     idealAngleAlt(k,1) = mod(subData.T.idealAngle(k) + angleShift,2*pi);
                 else
                     idealAngleAlt(k,1) = mod(subData.T.idealAngle(k) - angleShift,2*pi);
                 end
            elseif session == 3
                idealAngleAlt(k,1) = mod(subData.T.idealAngle(k) - angleShift,2*pi);
            end
                
        end
        
%   ---   TEST   replace actual angle with ideal angle for alternate state   
%     actualAngle = idealAngleAlt;
%     
% %     replace subset of trials with random responses
%     repPerRand = 0.95;
%     nreplaceRand = round(repPerRand*length(idealAngleAlt));
%     repIdx = randperm(length(idealAngleAlt));
%     randRepIdx = repIdx <= nreplaceRand;
%     actualAngle(randRepIdx) = rand(nreplaceRand,1)*2*pi-pi;
%     
% %     replace subset of trials with correct responses
%     repPerTarget = 0.05;
%     nreplaceTarget = round(repPerTarget*length(idealAngleAlt));
%     targetRepIdx = repIdx > nreplaceRand & repIdx <= nreplaceTarget + nreplaceRand;
%     actualAngle(targetRepIdx) = idealAngle(targetRepIdx);
%         
%  ----------------

    subData.T.actualAngle = actualAngle;
    subData.T.angleDiff = circ_dist(idealAngle,actualAngle);


    subData.T.idealAngleAlt = idealAngleAlt;
    subData.T.angleDiffAlt = abs((circ_dist((idealAngleAlt),(actualAngle))));
    subData.T.normError = abs(subData.T.angleDiff)./(abs(subData.T.angleDiff)+abs(subData.T.angleDiffAlt));


    % 

    if exist('groupData')
        groupData=[groupData; subData.T];
    else
        groupData=subData.T;
    end

end

% create a variable that corresponds to the previous mean:
guessPrevMean=nan(size(groupData.StimulusRep));




% 
if session == 1 || session == 3
    totStimRep=(groupData.Block-1).*8+groupData.StimulusRep;
elseif session == 2
    totStimRep=(groupData.Block-1).*6+groupData.StimulusRep;
end



% get angular distances between ideal and alternative for all trials
angleChange = (circ_dist(groupData.idealAngle,groupData.idealAngleAlt));




% calculate mean overall error across participants
SNuni = unique(groupData.SN);
for i = 1:length(SNuni)
    idx = groupData.SN==SNuni(i);
    meanErrorSs(i) = mean(abs(groupData.angleDiff(idx)));
end

% get median split for participants by mean absolute error...
medianErrorSs = meanErrorSs > median(meanErrorSs);
for i = 1:length(SNuni)
    idx = groupData.SN==SNuni(i);
    groupData.medianGroup(idx) = medianErrorSs(i);
end


% create selection array for transfer trials:
pStim=[0; groupData.StimulusRep(1:end-1)];
if session == 1 || session == 3
    sel = groupData.StimulusRep==1 & pStim==1 & groupData.Block>1;
elseif session == 2 && novelColors == false
    sel = groupData.transfer_trial & groupData.novel_color == 0;
elseif session == 2 && novelColors == true
    sel = groupData.transfer_trial & groupData.novel_color == 1 & groupData.ColorID ~= 13;
end

figure;
subplot(2, 1, 1)
hold on
histogram(groupData.angleDiff(sel), 50, 'Normalization','probability')
xlabel('Subject Error')
ylabel('Frequency')


subplot(2, 1, 2)


hold on
histogram(groupData.normError(sel), 50, 'Normalization','probability')
xlabel('norm error')
ylabel('Frequency')

figure;
subplot(2, 1, 1)
hold on
histogram(rad2deg(groupData.angleDiff(groupData.medianGroup == 1 & sel)), 50, 'Normalization','probability')
xlabel('Subject Error')
ylabel('Frequency')
title('> median');
xlim([-180,180])
xticks([-180:45:180]);


subplot(2, 1, 2)


hold on
histogram(rad2deg(groupData.angleDiff(groupData.medianGroup == 0 & sel)), 50, 'Normalization','probability')
xlabel('Subject Error')
ylabel('Frequency')
title('<= median');
xlim([-180,180])
xticks([-180:45:180]);



% loop through participants and compute angular error metric for transfer
% trials:
for i = 1:max(groupData.SN)
    subSel= sel& groupData.SN==i;
    %hist(groupData.angleDiff(subSel), 50)
    subErrNorm(i)=nanmean(abs(groupData.normError(subSel)));
    subErrAbsAngle(i)=nanmean(abs(groupData.angleDiff(subSel)));
%     prevErr(i)=nanmean(abs(guessPrevMean(subSel)));


    for j=1:nIt
    % create errors for a random model:
    rndErrNorm(i,j)= nanmean(abs(rand(sum(subSel), 1).*1-0));
    rndErrAbsAngle(i,j)= nanmean(abs(rand(sum(subSel), 1).*2.*pi-pi));
    end

end

% sort random guessing errors
sortedRndErrorNorm=sort(rndErrNorm)
maxErrorNorm = prctile(sortedRndErrorNorm',[97.5]);
minErrorNorm = prctile(sortedRndErrorNorm',[2.5]);

sortedRndErrorAbsAngle=sort(rndErrAbsAngle)
maxErrorAbsAngle = prctile(sortedRndErrorAbsAngle',[97.5]);
minErrorAbsAngle = prctile(sortedRndErrorAbsAngle',[2.5]);

% close all
figure;
hold on;
a=plot(sort(subErrNorm), 'r')
% b=plot(sort(prevErr), 'b')
plot(minErrorNorm, '--k')
plot(maxErrorNorm, '--k')

ff=legend([a], 'Participant error')
set(ff, 'location', 'northwest', 'box', 'off')
ylabel('normalized error')
xlabel('Participant rank')
ylim([0,1])


figure;
hold on;
a=plot(sort(subErrAbsAngle), 'r')
% b=plot(sort(prevErr), 'b')
plot(minErrorAbsAngle, '--k')
plot(maxErrorAbsAngle, '--k')

ff=legend([a], 'Participant error')
set(ff, 'location', 'northwest', 'box', 'off')
ylabel('absolute angular error')
xlabel('Participant rank')
ylim([0,pi]);
yticks([0,pi/4,pi/2,pi*0.75,pi]);
yticklabels([0,45,90,135,180]);
box on;


