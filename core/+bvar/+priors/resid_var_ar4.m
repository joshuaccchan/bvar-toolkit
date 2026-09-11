% bvar.priors.resid_var_ar4 - residual variances of univariate AR(4) models,
% used to set the Minnesota-prior scalings.
%
%   sig2 = bvar.priors.resid_var_ar4(Y0, Y)
%
%   Y0   : presample rows; the last 4 are prepended to Y
%   Y    : T x n estimation sample
%   sig2 : n-vector of AR(4) residual variances
%
% NEVER merge with bvar.priors.resid_var_allvars_ridge: that one regresses each
% variable on 4 lags of ALL variables with a 1e-4 ridge - numerically different
% sig2.
%
% See:
% Chan, J.C.C. (2021). Minnesota-Type Adaptive Hierarchical Priors for
% Large Bayesian VARs, International Journal of Forecasting, 37(3): 1212-1226

function sig2 = resid_var_ar4(Y0,Y)
[T,n] = size(Y);
sig2 = zeros(n,1);
tmpY = [Y0(end-4+1:end,:); Y];
for i=1:n
    Z = [ones(T,1) tmpY(4:end-1,i) tmpY(3:end-2,i) tmpY(2:end-3,i) tmpY(1:end-4,i)];
    tmpb = (Z'*Z)\(Z'*tmpY(5:end,i));
    sig2(i) = mean((tmpY(5:end,i)-Z*tmpb).^2);
end
end
