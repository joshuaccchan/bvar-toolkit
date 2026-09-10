% bvar.priors.resid_var_allvars_ridge - residual variances from regressing each
% variable on 4 lags of ALL variables (VAR(4)-style regressor set) with a 1e-4
% ridge on the normal equations; used to set the Minnesota-prior scalings in the
% hybrid TVP-VAR.
%
%   sig2 = bvar.priors.resid_var_allvars_ridge(Y0, Y)
%
%   Y0   : presample rows; the last 4 are prepended to Y
%   Y    : T x n estimation sample
%   sig2 : n-vector of residual variances
%
% NEVER merge with bvar.priors.resid_var_ar4: that one runs univariate AR(4)
% regressions with no ridge - numerically different sig2, hence different
% Minnesota scalings.
%
% Provenance and the legacy copies this stands in for: tests/variant_map.md.
% Equivalence: tests/unit/test_resid_var_allvars_ridge.m.
%
% See:
% Chan, J.C.C. (2023). Large Hybrid Time-Varying Parameter VARs, Journal of
% Business and Economic Statistics, 41(3): 890-905

function sig2 = resid_var_allvars_ridge(Y0,Y)
[T,n] = size(Y);
sig2 = zeros(n,1);
tmpY = [Y0(end-4+1:end,:); Y];
for i=1:n
    Z = [ones(T,1) tmpY(4:end-1,:) tmpY(3:end-2,:) tmpY(2:end-3,:) tmpY(1:end-4,:)];
    tmpb = (Z'*Z+1e-4*eye(size(Z,2)))\(Z'*tmpY(5:end,i));
    sig2(i) = mean((tmpY(5:end,i)-Z*tmpb).^2);
end
end
