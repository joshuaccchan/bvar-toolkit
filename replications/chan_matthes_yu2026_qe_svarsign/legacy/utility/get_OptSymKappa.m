% This function obtains the optimal shrinkage hyperparameter values and the
% associated marginal likelihood value under the symmetric prior
%
% See:
% Chan, J.C.C. (2021). Asymmetric conjugate priors for large Bayesian VARs,
% Quantitative Economics, forthcoming.
%
% Input: idx_ns - index for nonstationary variables

function [ml_opt,kappa_opt] = get_OptSymKappa(Y0,Y,Z,p,type,idx_ns)
if nargin == 5
    idx_ns = [];
end
kappa3 = 1; kappa4 = 100;
n = size(Y,2);
sig2 = get_resid_var(Y0,Y);
if strcmp(type,'stru')
    f = @(k1) -ml_VAR_ACP(p,Y,Z,prior_ACP_stru(n,p,[k1,k1,kappa3,kappa4],sig2,idx_ns));
elseif strcmp(type,'redu')
    f = @(k1) -ml_VAR_ACP(p,Y,Z,prior_ACP_redu(n,p,[k1,k1,kappa3,kappa4],sig2,idx_ns));
end
[kappa1,nml] = fminbnd(f,0,1);
ml_opt = -nml;
kappa_opt = [kappa1,kappa1,kappa3,kappa4];
end