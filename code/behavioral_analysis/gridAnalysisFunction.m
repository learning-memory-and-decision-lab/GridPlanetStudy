function [T, store] = gridAnalysisFunction(subjNum,ThreeRunSubj)
%make a table for the behavioral data for fmri (session 2) of gridPlanet task
%will output the behavioral table for the data + the store struct for fmri
%timing
datadir = ['/oscar/data/mnassar/lyu21/GridPlanet/behavioralData'];

if sum(subjNum == ThreeRunSubj) > 0
   load(sprintf('%s/gridplanet_Scanning_subj_%d_block_3.mat', datadir,subjNum))
else
   load(sprintf('%s/gridplanet_Scanning_subj_%d_block_4.mat', datadir,subjNum))
end

% if this is some particular subjects where the task was restarted,
% separately load their colours2angles var from earlier session
if subjNum == 1065
    load(sprintf('%s/gridplanet_Scanning_subj_%d_block_3.mat', datadir,subjNum),'colours2Angles');
end

numTrials = length(Trial_accuracy);
numColors = length(colours2Angles);

    numaddColors = size(addcolours2Angles,1);


 %give colours a name
for t=1:length(Trial_colours2Angles)
    
    % if this is some particular subjects where the task was restarted,
    % separately load their colours2angles var from earlier session
    if (subjNum == 1065) && t <= 252
        load(sprintf('%s/gridplanet_Scanning_subj_%d_block_3.mat', datadir,subjNum),'colours2Angles');
    elseif (subjNum == 1065) && t > 252
        load(sprintf('%s/gridplanet_Scanning_subj_%d_block_4.mat', datadir,subjNum),'colours2Angles');
    end
    if (subjNum == 1074) && t <= 168
        load(sprintf('%s/gridplanet_Scanning_subj_%d_block_2.mat', datadir,subjNum),'colours2Angles');
    elseif (subjNum == 1074) && t > 168
        load(sprintf('%s/gridplanet_Scanning_subj_%d_block_4.mat', datadir,subjNum),'colours2Angles');
    end
    
    if Trial_colours2Angles(1:3,t)' == colours2Angles(1,1:3)
        Trial_colours2Angles(6,t) = 1;
    elseif Trial_colours2Angles(1:3,t)' == colours2Angles(2,1:3)
        Trial_colours2Angles(6,t) = 2;
    elseif Trial_colours2Angles(1:3,t)' == colours2Angles(3,1:3)
         Trial_colours2Angles(6,t) = 3;
    elseif Trial_colours2Angles(1:3,t)' == colours2Angles(4,1:3)
         Trial_colours2Angles(6,t) = 4;
    elseif Trial_colours2Angles(1:3,t)' == colours2Angles(5,1:3)
         Trial_colours2Angles(6,t) = 5;
    end

  
 
        if Trial_colours2Angles(1:3,t)' == addcolours2Angles(1,1:3)
            Trial_colours2Angles(6,t) = 6; %set1
        elseif Trial_colours2Angles(1:3,t)' == addcolours2Angles(2,1:3)
            Trial_colours2Angles(6,t) = 7;%set1
        elseif Trial_colours2Angles(1:3,t)' == addcolours2Angles(3,1:3)
            Trial_colours2Angles(6,t) = 8;%set1
        elseif Trial_colours2Angles(1:3,t)' == addcolours2Angles(4,1:3)
            Trial_colours2Angles(6,t) = 9;%set1
        elseif Trial_colours2Angles(1:3,t)' == addcolours2Angles(5,1:3)
            Trial_colours2Angles(6,t) = 10; %set 2
        elseif Trial_colours2Angles(1:3,t)' == addcolours2Angles(6,1:3)
            Trial_colours2Angles(6,t) = 11; %set 2
        elseif Trial_colours2Angles(1:3,t)' == addcolours2Angles(7,1:3)
            Trial_colours2Angles(6,t) = 12; %set 2
        elseif Trial_colours2Angles(1:3,t)' == addcolours2Angles(8,1:3)
            Trial_colours2Angles(6,t) = 13; %set 2
            
        end
  
  
end

if sum(subjNum == ThreeRunSubj) > 0
    blockID = [ones(numLocations*numStimReps,1); ones(numLocations*numStimReps,1)*2; ones(numLocations*numStimReps,1)*3; ...
    ones(numLocations*numStimReps,1)*4; ones(numLocations*numStimReps,1)*5; ...
    ones(numLocations*numStimReps,1)*6];
    MRIblockID = [ones(numLocations*numStimReps,1); ones(numLocations*numStimReps,1)*1; ones(numLocations*numStimReps,1)*2; ...
    ones(numLocations*numStimReps,1)*2; ones(numLocations*numStimReps,1)*3; ...
    ones(numLocations*numStimReps,1)*3];

else
    blockID = [ones(numLocations*numStimReps,1); ones(numLocations*numStimReps,1)*2; ones(numLocations*numStimReps,1)*3; ...
        ones(numLocations*numStimReps,1)*4; ones(numLocations*numStimReps,1)*5; ...
        ones(numLocations*numStimReps,1)*6; ones(numLocations*numStimReps,1)*7; ones(numLocations*numStimReps,1)*8];
    MRIblockID = [ones(numLocations*numStimReps,1); ones(numLocations*numStimReps,1)*1; ones(numLocations*numStimReps,1)*2; ...
        ones(numLocations*numStimReps,1)*2; ones(numLocations*numStimReps,1)*3; ...
        ones(numLocations*numStimReps,1)*3; ones(numLocations*numStimReps,1)*4; ones(numLocations*numStimReps,1)*4];
end

%colorBlockID = [ones(120,1); ones(120,1)*2; ones(120,1)*3; ones(120,1)*4];




     T = table(blockID, MRIblockID, Trial_stimRep', Trial_WaitSim', Trial_colours2Angles(6,:)', ...
   Trial_colours2Angles(4,:)', Trial_colours2Angles(5,:)', Trial_startingshieldx', Trial_startingshieldy', ...
    Trial_endshieldx', Trial_endshieldy', destx', desty', Trial_responseAngle', Trial_accuracy', RTs', 'VariableNames', ...
    {'Block', 'MRIBlock','StimulusRep', 'WaitSimTime', 'ColorID', ...
      'MeanX', 'MeanY', 'ShieldStartX', 'ShieldStartY', 'RespX', 'RespY', ...
      'OutcomeX', 'OutcomeY','respAngle', 'Accuracy', 'RT'});
  filename2 = sprintf('sortedData_subject%d_session2', subjNum);
   % save(filename2, 'T')


end



