% This script obtains the optimal kappa1 and kappa2 by maximizing the
% marginal likelihood, assuming kappa1 = kappa2
%
% See:
% Chan, J.C.C. (2019). Asymmetric conjugate priors for large Bayesian VARs,
% CAMA Working Papers 51/2019

function [ml_opt,kappa_opt] = get_OptSymKappa(Y0,Y,Z,p)
kappa3 = 1; kappa4 = 100;
f = @(k1) -ml_VAR_ACP(p,Y,Y0,Z,[k1,k1,kappa3,kappa4]);
[kappa1,nml] = fminbnd(f,0,1);
ml_opt = -nml;
kappa_opt = [kappa1,kappa1,kappa3,kappa4];
end