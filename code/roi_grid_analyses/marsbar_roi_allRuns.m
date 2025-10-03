function [b,roi_file]=marsbar_roi_allRuns(subj,glmPath,mask)
% this function gets the percent signal change from ROIs specified in a
% path for specified SPM.mat files. 
% code is taken from marsbar help site
% (http://marsbar.sourceforge.net/faq.html)
% modified by AV, 2018

% output structure
 outputdir = sprintf('/gpfs/data/mnassar/lyu21/GridPlanet/allrunsGLMs/%s/%d/',glmPath,subj);

for i = 1:length(subj)
    % specify directory
    spm_name = sprintf('/gpfs/data/mnassar/lyu21/GridPlanet/allrunsGLMs/%s/%d/SPM.mat',glmPath, subj);
    %    spm_name = sprintf('/gpfs/data/mnassar/lyu21/GridPlanet/allrunsGLMs/stimMod_dropErrors_last2runs/%d/SPM.mat', subj);

    if strcmp(mask,'ppc')
     roi_file = '/gpfs/data/mnassar/lyu21/GridPlanet/Templates/coggrid_allruns_ppc_fullbatch_sphere_gm_roi.mat';
    elseif strcmp(mask,'mpfc_allgrids')
     roi_file = '/gpfs/data/mnassar/lyu21/GridPlanet/Templates/allgrid_mpfc_fullbatch_sphere_gm_roi.mat';
    elseif strcmp(mask,'pcc_allgrids')
     roi_file = '/gpfs/data/mnassar/lyu21/GridPlanet/Templates/allgrid_pcc_fullbatch_sphere_gm_roi.mat';
    elseif strcmp(mask,'ltpj_allgrids')
     roi_file = '/gpfs/data/mnassar/lyu21/GridPlanet/Templates/allgrid_ltpj_fullbatch_sphere_gm_roi.mat';
    elseif strcmp(mask,'rtpj_allgrids')
     roi_file = '/gpfs/data/mnassar/lyu21/GridPlanet/Templates/allgrid_rtpj_fullbatch_sphere_gm_roi.mat';
    elseif strcmp(mask,'rlpfc_allgrids')
     roi_file = '/gpfs/data/mnassar/lyu21/GridPlanet/Templates/allgrid_rlpfc_fullbatch_sphere_gm_roi.mat';
    elseif strcmp(mask,'constantinescu_tpj')
     roi_file = '/gpfs/data/mnassar/lyu21/GridPlanet/Templates/constantinescu_tpj_gm_roi.mat';
    elseif strcmp(mask,'constantinescu_vmpfc')
     roi_file = '/gpfs/data/mnassar/lyu21/GridPlanet/Templates/constantinescu_vmpfc_gm_roi.mat';
    elseif strcmp(mask,'constantinescu_pcc')
     roi_file = '/gpfs/data/mnassar/lyu21/GridPlanet/Templates/constantinescu_pcc_gm_roi.mat';
    elseif strcmp(mask,'constantinescu_ppc')
     roi_file = '/gpfs/data/mnassar/lyu21/GridPlanet/Templates/constantinescu_ppc_gm_roi.mat';
    elseif strcmp(mask,'erc_left')
     roi_file = '/gpfs/data/mnassar/lyu21/GridPlanet/Templates/pmEC_PHCpref_left_MNI_roi.mat';
    elseif strcmp(mask,'erc_right')
     roi_file = '/gpfs/data/mnassar/lyu21/GridPlanet/Templates/pmEC_PHCpref_right_MNI_roi.mat';
    elseif strcmp(mask,'erc_right_park')
     roi_file = '/gpfs/data/mnassar/lyu21/GridPlanet/Templates/park_ECr_roi.mat';
    elseif strcmp(mask,'erc_left_park')
     roi_file = '/gpfs/data/mnassar/lyu21/GridPlanet/Templates/park_ECl_roi.mat';
    end

    % Make marsbar design object
    D  = mardo(spm_name);
    % Make marsbar ROI object
    R  = maroi(roi_file);
    % Fetch data into marsbar data object
    Y  = get_marsy(R, D, 'mean');
    % Get contrasts from original design
  %  xCon = get_contrasts(D);
    % Estimate design on ROI data
   % r = residuals(D);
    E = estimate(D, Y);
    % Put contrasts from original design back into design object
   % E = set_contrasts(E, xCon);
    % get design betas
    b = betas(E);
    % get stats and stuff for all contrasts into statistics structure
  %  marsS = compute_contrasts(E, 1:length(xCon));
     filename = [outputdir sprintf('betas_dropErrors_%s', mask)];
     save(filename, 'b');

%     end
end