% bvar.priors.acp_opt_kappa - choose the two shrinkage hyperparameters of the
% asymmetric conjugate prior by maximizing the closed-form marginal likelihood.
%
%   [ml_opt,kappa_opt] = bvar.priors.acp_opt_kappa(Y0, Y, Z, p, k0, type)
%   [ml_opt,kappa_opt] = bvar.priors.acp_opt_kappa(Y0, Y, Z, p, k0, type, idx_ns)
%   [...]              = bvar.priors.acp_opt_kappa(..., 'symmetric', true)
%
%   Y0, Y  : presample rows and the T x n estimation sample
%   Z      : T x (n p + 1) lag matrix, intercept first (bvar.util.build_lags)
%   p      : lag length
%   k0     : starting values [kappa1, kappa2]; ignored when 'symmetric' is true
%   type   : 'redu' or 'stru', selecting bvar.priors.acp_redu or acp_stru
%   idx_ns : indices of variables entering in levels (default none)
%   'symmetric' : default false, which optimizes log kappa1 and log kappa2 with
%               fminsearch, so the arguments stay positive without a constrained
%               solver; true imposes kappa1 = kappa2 and uses fminbnd on (0,1)
%   ml_opt    : the maximized log marginal likelihood, bvar.ml.acp at kappa_opt
%   kappa_opt : the full 4-vector [kappa1, kappa2, 1, 100] - kappa3 and kappa4
%               are held at those values, as in the paper's application
%
% Because bvar.ml.acp is available in closed form, this is a two-parameter
% optimization over a smooth objective, with no repeated estimation of the VAR.
%
% Provenance and the legacy copies this stands in for: tests/variant_map.md.
% Equivalence: tests/unit/test_acp_equivalence.m.
%
% See:
% Chan, J.C.C. (2022). Asymmetric Conjugate Priors for Large Bayesian VARs,
% Quantitative Economics, 13(3): 1145-1169.

function [ml_opt,kappa_opt] = acp_opt_kappa(Y0,Y,Z,p,k0,type,idx_ns,varargin)
if nargin < 7 || isempty(idx_ns)
    idx_ns = [];
end
symmetric = false;
for iv = 1:2:numel(varargin)
    switch lower(varargin{iv})
        case 'symmetric', symmetric = varargin{iv+1};
        otherwise, error('bvar:priors:acp_opt_kappa:badOption', ...
                'unknown option ''%s''', varargin{iv});
    end
end
assert(any(strcmp(type,{'redu','stru'})), 'bvar:priors:acp_opt_kappa:badType', ...
    'type must be ''redu'' or ''stru'', got ''%s''', type);

kappa3 = 1; kappa4 = 100;
n = size(Y,2);
sig2 = bvar.priors.resid_var_ar4(Y0,Y);

if symmetric
    if strcmp(type,'stru')
        f = @(k1) -bvar.ml.acp(p,Y,Z,bvar.priors.acp_stru(n,p,[k1,k1,kappa3,kappa4],sig2,idx_ns));
    else
        f = @(k1) -bvar.ml.acp(p,Y,Z,bvar.priors.acp_redu(n,p,[k1,k1,kappa3,kappa4],sig2,idx_ns));
    end
    [kappa1,nml] = fminbnd(f,0,1);
    ml_opt = -nml;
    kappa_opt = [kappa1,kappa1,kappa3,kappa4];
else
    if strcmp(type,'stru')
        f = @(k) -bvar.ml.acp(p,Y,Z,bvar.priors.acp_stru(n,p,[exp(k(1)),exp(k(2)),kappa3,kappa4],sig2,idx_ns));
    else
        f = @(k) -bvar.ml.acp(p,Y,Z,bvar.priors.acp_redu(n,p,[exp(k(1)),exp(k(2)),kappa3,kappa4],sig2,idx_ns));
    end
    [k_opt,nml] = fminsearch(f,log(k0));
    ml_opt = -nml;
    kappa_opt = [exp(k_opt(1)),exp(k_opt(2)),kappa3,kappa4];
end
end
