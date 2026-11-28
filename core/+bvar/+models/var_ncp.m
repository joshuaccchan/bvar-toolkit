% bvar.models.var_ncp - independent posterior draws of a VAR under the natural conjugate prior,
% with the shrinkage chosen by the marginal likelihood. Written for this toolkit.
%
%   res = bvar.models.var_ncp(Y0, Y, p)
%   res = bvar.models.var_ncp(Y0, Y, p, 'lev', [3 4], 'c1', .04, 'nsim', 0)
%
%   Y0 : presample, at least max(p, 4) rows, same columns as Y
%   Y  : T x n sample
%   p  : lag length
%
%   'lev'  : the variables in levels, whose own first lag has prior mean one; default none
%   'c1'   : shrinkage of the lags; default the value that maximizes the marginal likelihood
%            over [1e-3, 10]
%   'c3'   : prior variance of the intercepts, default 100
%   'nsim' : number of draws, default 10000; 0 returns no draws
%   'omit' : indices of rows of Y left out of the likelihood, such as the months of the pandemic;
%            they still serve as lags of the rows after them, and the AR(4) fits for s2 use
%            every row; default none
%
%   res.c1, res.lml : the shrinkage and the log marginal likelihood at it
%   res.s2  : n x 1 residual variances of AR(4) fits (bvar.priors.resid_var_ar4)
%   res.VA, res.B0 : k x 1 prior variances and k x n prior means of the coefficients on
%             x_t = (1, y_{t-1}', ..., y_{t-p}')', k = 1 + n*p
%   res.A   : nsim x k*n coefficients, each row a k x n matrix A(:)'
%   res.Sig : nsim x n x n error covariance matrices
%
% The prior is vec(A) | Sig ~ N(vec(B0), Sig kron diag(VA)) and Sig ~ IW(n + 2, diag(s2)), with
% prior variance c3 for the intercepts and c1/(l^2 s2(r)) for lag l of variable r. The draws are
% independent. res.A and res.Sig have the layout of bvar.forecast.predictive.
%
% See:
% Chan, J.C.C. (2027). Bayesian Macroeconometrics: Methods and Applications. Chapman &
% Hall/CRC.

function res = var_ncp(Y0, Y, p, varargin)
lev = []; c1 = []; c3 = 100; nsim = 10000; omit = [];
iv = 1;
while iv <= numel(varargin)
    if iv == numel(varargin)
        error('bvar:models:var_ncp:badOption', 'option ''%s'' has no value', ...
            char(string(varargin{iv})));
    end
    switch lower(char(string(varargin{iv})))
        case 'lev',  lev = varargin{iv+1};
        case 'c1',   c1 = varargin{iv+1};
        case 'c3',   c3 = varargin{iv+1};
        case 'nsim', nsim = varargin{iv+1};
        case 'omit', omit = varargin{iv+1};
        otherwise, error('bvar:models:var_ncp:badOption', 'unknown option ''%s''', ...
                char(string(varargin{iv})));
    end
    iv = iv + 2;
end

[T, n] = size(Y);
if size(Y0, 2) ~= n
    error('bvar:models:var_ncp:badData', 'Y0 has %d columns and Y has %d', size(Y0, 2), n);
end
if size(Y0, 1) < max(p, 4)
    error('bvar:models:var_ncp:badData', ...
        'Y0 has %d rows; the lags and the AR(4) fits need at least %d', size(Y0, 1), max(p, 4));
end
if ~all(isfinite([Y0(:); Y(:)]))
    error('bvar:models:var_ncp:badData', 'Y0 and Y contain NaN or Inf');
end
if ~all(ismember(lev, 1:n)) || (~isempty(c1) && ~(isscalar(c1) && c1 > 0)) ...
        || ~(isscalar(c3) && c3 > 0) || ~(isscalar(nsim) && nsim >= 0 && nsim == fix(nsim))
    error('bvar:models:var_ncp:badOption', ['lev must index columns of Y, c1 and c3 must be ' ...
        'positive and nsim a nonnegative integer']);
end
if ~isempty(omit) && ~(isnumeric(omit) && all(ismember(omit(:), 1:T)) && numel(unique(omit)) < T)
    error('bvar:models:var_ncp:badOption', 'omit must index rows of Y and leave some');
end

[~, Z] = bvar.util.build_lags([Y0(end-p+1:end, :); Y], p);
k = size(Z, 2);
s2 = bvar.priors.resid_var_ar4(Y0, Y);
keep = true(T, 1); keep(omit) = false;
Y = Y(keep, :); Z = Z(keep, :); T = size(Y, 1);
if isempty(c1)
    [lc, f] = fminbnd(@(lc) -ncp_lml(Y, Z, p, lev, exp(lc), c3, s2), log(1e-3), log(10));
    c1 = exp(lc); lml = -f;
else
    lml = ncp_lml(Y, Z, p, lev, c1, c3, s2);
end
[VA, B0] = ncp_prior(p, n, lev, c1, c3, s2);

A = zeros(nsim, k*n); Sig = zeros(nsim, n, n);
if nsim > 0
    [~, Bhat, K, Shat] = bvar.ml.ncp(Y, Z, VA, n + 2, diag(s2), B0);
    C = chol(K, 'lower');
    for r = 1:nsim
        Sigr = iwishrnd(Shat, n + 2 + T);
        Ar = Bhat + (C'\randn(k, n))*chol(Sigr, 'lower')';
        A(r, :) = Ar(:)';
        Sig(r, :, :) = Sigr;
    end
end
res = struct('c1', c1, 'lml', lml, 's2', s2, 'VA', VA, 'B0', B0, 'A', A, 'Sig', Sig);
end

function [VA, B0] = ncp_prior(p, n, lev, c1, c3, s2)
% the prior variances VA and means B0 of the coefficients on x_t: a mean of one on the own first
% lag of the variables lev, zero elsewhere
VA = [c3; kron(1./(1:p)'.^2, c1./s2)];
B0 = zeros(1 + n*p, n);
B0(sub2ind(size(B0), 1 + lev, lev)) = 1;
end

function lml = ncp_lml(Y, Z, p, lev, c1, c3, s2)
% the log marginal likelihood at c1
n = size(Y, 2);
[VA, B0] = ncp_prior(p, n, lev, c1, c3, s2);
lml = bvar.ml.ncp(Y, Z, VA, n + 2, diag(s2), B0);
end
