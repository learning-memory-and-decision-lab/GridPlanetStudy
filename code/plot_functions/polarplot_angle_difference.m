
% figure for angle differences between conditions
clear all

% load the betas
% load('/oscar/data/mnassar/lyu21/GridPlanet/roi_analysis/phaseConsistency_state_allruns/pcc_gm_allgrids.mat','slist','sin_cos_betas')

numFolds = 6;
within = false;


if ~within
    load('/oscar/data/mnassar/lyu21/GridPlanet/roi_analysis/phaseConsistency_state_allruns/ppc_gm_coggrid.mat','slist','sin_cos_betas')

    for i = 1:length(slist)
        phi0(i,1) = mod(atan2(sin_cos_betas(i,1),sin_cos_betas(i,2)),pi*2)./numFolds;
        phi1(i,1) = mod(atan2(sin_cos_betas(i,3),sin_cos_betas(i,4)),pi*2)./numFolds;
    end
    dist = abs(phi0-phi1);

else
    load('/oscar/data/mnassar/lyu21/GridPlanet/roi_analysis/phaseConsistency_withinstate/ppc_gm_withinstate.mat','slist','sin_cos_betas')

    sin_cos_betas(sin_cos_betas==0)=NaN;
    for i = 1:length(slist)
        for b = 1:4
            phi0_b(i,b) = mod(atan2(sin_cos_betas(i,b,1),sin_cos_betas(i,b,2)),pi*2);
            phi1_b(i,b) = mod(atan2(sin_cos_betas(i,b,3),sin_cos_betas(i,b,4)),pi*2);
        end

        % calculate difference WITHIN states
        if slist(i) == 1066 || slist(i) == 1080
            nrun = 3;
        else
            nrun = 4;
        end
        pairs = nchoosek(1:nrun,nrun-1);
        diff0_p = nan;
        diff1_p = nan;
        for p = 1:size(pairs,1)
            runOut = find(sum(1:nrun~=pairs(p,:)')==nrun-1);
            meanIn0 = circ_mean(phi0_b(i,pairs(p,:)),[],2);
            meanIn1 = circ_mean(phi1_b(i,pairs(p,:)),[],2);
            diff0_p(p,1) = circ_dist(phi0_b(i,runOut),meanIn0);
            diff1_p(p,1) = circ_dist(phi1_b(i,runOut),meanIn1);
        end

        % get mean
        idx0 = ~isnan(diff0_p);
        idx1 = ~isnan(diff1_p);

        diff0 = circ_mean(diff0_p(idx0));
        diff1 = circ_mean(diff1_p(idx1));

        dist(i,1) = circ_mean([diff0;diff1])./numFolds;


    end

end


figure;polarhistogram(6*dist,20);
thetaticks(0:90:360)
thetaticklabels([0,15,30,45])

