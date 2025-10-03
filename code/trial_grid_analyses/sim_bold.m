function [convresp,r] = sim_bold(TR,startrespTR,endrespTR,values,maxduration)
% simulate bold signal using inputted trial times in TRs (trialTimesTR) and
% some values in a vector of the same length (values), maxduration
% specifies end of simulated vector in TRs

% specify microtime
dt = 0.05;

% make vector for plopping in values
veclength = (maxduration*TR)/dt;
R = zeros(round(veclength),1);

% %           get basis fxn
xBF.dt = dt;
xBF.name = 'hrf';
bf = spm_get_bf(xBF);
%             specify timings and amplitudes
idxStart = round((startrespTR*TR)/dt);
idxEnd = round((endrespTR*TR)/dt);

r = zeros(size(R,1),1);
r(idxStart:idxEnd) = values;

% bf = spm_hrf(TR);
% convresp = conv(bf,r);
% convresp = convresp(1:maxduration);

% %             specify convolution
U.u = r;
U.name = {'reg'};
%             convolve rt with hrf
convresp = spm_Volterra(U,bf.bf);

% resample time series
convresp = convresp(1:(TR/dt):veclength);
