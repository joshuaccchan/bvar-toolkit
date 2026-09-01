% This script obtains the optimal kappa1 and kappa2 by maximizing the
% marginal likelihood
%
% See:
% Chan, J.C.C. (2019). Asymmetric conjugate priors for large Bayesian VARs,
% CAMA Working Papers 51/2019

function [ml_opt,kappa_opt] = get_OptKappa(Y0,Y,Z,p,k0)
kappa3 = 1; kappa4 = 100;
f = @(k) -ml_VAR_ACP(p,Y,Y0,Z,[exp(k(1)),exp(k(2)),kappa3,kappa4]);
[k_opt,nml] = fminsearch(f,log(k0));
ml_opt = -nml;
kappa_opt = [exp(k_opt(1)),exp(k_opt(2)),kappa3,kappa4];
end