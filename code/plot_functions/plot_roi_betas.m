function plot_roi_betas
% this is a function for plotting beta coefficients from CV analysis of
% state-state grid codes

% folder to get coefficients
datadir = '/oscar/data/mnassar/lyu21/GridPlanet/roi_analysis/phaseConsistency_state_allruns/';

% roi names
roi = {'erc_left','erc_right','ltpj_gm_allgrids','rtpj_gm_allgrids',...
    'ppc_gm_coggrid','pcc_gm_allgrids','mpfc_gm_allgrids'};

% load coefficients
for i = 1:length(roi)
    roi(i).data = load(fullfile(datadir,sprintf('%s.mat',roi{i})));
    beta(:,i) = roi(i).data.gmean_b(:,1);
end
