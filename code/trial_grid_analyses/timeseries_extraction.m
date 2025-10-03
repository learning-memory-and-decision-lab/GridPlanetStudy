function [y,design_obj] = timeseries_extraction(spmSubjPath,roiPath)

% % base directory
% basedir = '/gpfs/data/mnassar/lyu21/GridPlanet/';
% 
% % model name
% modelname = 'stimMod';
% 
% % Path for the SPM.mat file
% spmSubjPath = fullfile(basedir, int2str(subjectnum), ...
%     modelname, 'SPM.mat');

% 
% % Load ROI marsbar object
% roiFPath = fullfile('/Users/badrelab/Documents/RestaurantTask_fMRI/roi/', 'vmpfc_masked_bartra03D_roi.mat');

% Loads an object named 'roi' created in the marsbar toolbox GUI
load(roiPath);

% Make marsbar design object from the SPM.mat containing bold info
design_obj = mardo(spmSubjPath);

% load SPM mat file
load(spmSubjPath);

% Get data
data_obj = get_marsy(roi,design_obj,'mean');

% flags
flags.sessions=[1];



% apply filter to data
Y = apply_filter(design_obj, data_obj);



% get summary time course(s)
y = summary_data(Y); 

% % first TR to use 
% sTR=0;
% 
% % mean-center each run
% y1 = zscore(y(sTR+1:SPM.nscan(1)));
% y2 = zscore(y(SPM.nscan(1)+sTR+1:sum(SPM.nscan(1:2))));
% y3 = zscore(y(sum(SPM.nscan(1:2))+sTR+1:sum(SPM.nscan(1:3))));
% 
% % put data into a cell
% Y={};
% Y{1}=y1;
% Y{2}=y2;
% Y{3}=y3;
% % concatinate data
% y = [y1;y2;y3];

end
