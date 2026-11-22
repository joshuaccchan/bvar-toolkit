% bvar.models.var_soe - posterior draws of a VAR for a small open economy, in which the domestic
% variables do not enter the foreign equations (block exogeneity). Each block has a natural
% conjugate prior, and the draws come block by block. Written for this toolkit.
%
%   res = bvar.models.var_soe(Y0, Y, p, nf)
%   res = bvar.models.var_soe(Y0, Y, p, nf, 'lev', [3 4 7], 'c1', .04, 'nsim', 0)
%
%   Y0 : presample, at least max(p, 4) rows, same columns as Y
%   Y  : T x n sample, the nf foreign variables first
%   p  : lag length
%   nf : number of foreign variables
%
%   'lev'  : the variables in levels, whose own first lag has prior mean one; default none
%   'c1'   : shrinkage of the lags; default the value that maximizes the marginal likelihood
%            over [1e-3, 10]
%   'c3'   : prior variance of the intercepts, default 100
%   'nsim' : number of draws, default 10000; 0 returns no draws
%
%   res.c1, res.lml : the shrinkage and the log marginal likelihood at it
%   res.s2  : n x 1 residual variances of AR(4) fits (bvar.priors.resid_var_ar4)
%   res.VA, res.B0 : k x 1 prior variances and k x n prior means of the coefficients on
%             x_t = (1, y_{t-1}', ..., y_{t-p}')', k = 1 + n*p
%   res.A   : nsim x k*n reduced-form coefficients, each row a k x n matrix A(:)'
%   res.Sig : nsim x n x n reduced-form error covariance matrices
%
% The prior variance is c3 for the intercepts and c1/(l^2 s2(r)) for lag l of variable r. The
% domestic block also has the current foreign variables, with prior variances 1/s2(r). The error
% covariance matrices are IW(nf + 2, diag(s2(1:nf))) for the foreign block and IW(n + 2,
% diag(s2(nf+1:n))) for the domestic block, which gives the reduced-form error covariance matrix
% the prior IW(n + 2, diag(s2)) of an unrestricted VAR (Chan, 2022). The draws are independent.
% res.A and res.Sig have the layout of bvar.forecast.predictive.
%
% See:
% Cushman, D.O. and Zha, T. (1997). Identifying Monetary Policy in a Small Open Economy under
% Flexible Exchange Rates, Journal of Monetary Economics, 39(3): 433-448.
% Zha, T. (1999). Block Recursion and Structural Vector Autoregressions, Journal of Econometrics,
% 90(2): 291-316.
% Chan, J.C.C. (2022). Asymmetric Conjugate Priors for Large Bayesian VARs, Quantitative
% Economics, 13(3): 1145-1169.

function res = var_soe(Y0, Y, p, nf, varargin)
lev = []; c1 = []; c3 = 100; nsim = 10000;
iv = 1;
while iv <= numel(varargin)
    if iv == numel(varargin)
        error('bvar:models:var_soe:badOption', 'option ''%s'' has no value', ...
            char(string(varargin{iv})));
    end
    switch lower(char(string(varargin{iv})))
        case 'lev',  lev = varargin{iv+1};
        case 'c1',   c1 = varargin{iv+1};
        case 'c3',   c3 = varargin{iv+1};
        case 'nsim', nsim = varargin{iv+1};
        otherwise, error('bvar:models:var_soe:badOption', 'unknown option ''%s''', ...
                char(string(varargin{iv})));
    end
    iv = iv + 2;
end

[T, n] = size(Y);
if ~(isscalar(nf) && nf >= 1 && nf < n && nf == fix(nf))
    error('bvar:models:var_soe:badData', 'nf must be an integer from 1 to %d', n - 1);
end
if size(Y0, 2) ~= n
    error('bvar:models:var_soe:badData', 'Y0 has %d columns and Y has %d', size(Y0, 2), n);
end
if size(Y0, 1) < max(p, 4)
    error('bvar:models:var_soe:badData', ...
        'Y0 has %d rows; the lags and the AR(4) fits need at least %d', size(Y0, 1), max(p, 4));
end
if ~all(isfinite([Y0(:); Y(:)]))
    error('bvar:models:var_soe:badData', 'Y0 and Y contain NaN or Inf');
end
if ~all(ismember(lev, 1:n)) || (~isempty(c1) && ~(isscalar(c1) && c1 > 0)) ...
        || ~(isscalar(c3) && c3 > 0) || ~(isscalar(nsim) && nsim >= 0 && nsim == fix(nsim))
    error('bvar:models:var_soe:badOption', ['lev must index columns of Y, c1 and c3 must be ' ...
        'positive and nsim a nonnegative integer']);
end

nd = n - nf; fo = 1:nf; dm = nf+1:n;   % the foreign and the domestic variables
[~, Z] = bvar.util.build_lags([Y0(end-p+1:end, :); Y], p);
k = size(Z, 2);
iff = [1, reshape(1 + fo' + (0:p-1)*n, 1, [])];   % the intercept and the foreign lags in x_t
s2 = bvar.priors.resid_var_ar4(Y0, Y);
if isempty(c1)
    [lc, f] = fminbnd(@(lc) -soe_lml(Y, Z, nf, lev, exp(lc), c3, s2), log(1e-3), log(10));
    c1 = exp(lc); lml = -f;
else
    lml = soe_lml(Y, Z, nf, lev, c1, c3, s2);
end
[VA, B0] = soe_prior(p, n, lev, c1, c3, s2);

A = zeros(nsim, k*n); Sig = zeros(nsim, n, n);
if nsim > 0
    [~, Bf, Kf, Sf] = bvar.ml.ncp(Y(:, fo), Z(:, iff), VA(iff), nf + 2, diag(s2(fo)), B0(iff, fo));
    [~, Bd, Kd, Sd] = bvar.ml.ncp(Y(:, dm), [Z, Y(:, fo)], [VA; 1./s2(fo)], n + 2, ...
        diag(s2(dm)), [B0(:, dm); zeros(nf, nd)]);
    Cf = chol(Kf, 'lower'); Cd = chol(Kd, 'lower');
    for r = 1:nsim
        Sfr = iwishrnd(Sf, nf + 2 + T);                  % the foreign block
        Bfr = Bf + (Cf'\randn(numel(iff), nf))*chol(Sfr, 'lower')';
        Omr = iwishrnd(Sd, n + 2 + T);                   % the domestic block, independent of it
        Bdr = Bd + (Cd'\randn(k + nf, nd))*chol(Omr, 'lower')';
        [Ar, Sigr] = bvar.util.soe_reduced_form(Bfr, Bdr, Sfr, Omr);
        A(r, :) = Ar(:)';
        Sig(r, :, :) = Sigr;
    end
end
res = struct('c1', c1, 'lml', lml, 's2', s2, 'VA', VA, 'B0', B0, 'A', A, 'Sig', Sig);
end

function [VA, B0] = soe_prior(p, n, lev, c1, c3, s2)
% the prior variances VA and means B0 of the coefficients on x_t: a mean of one on the own first
% lag of the variables lev, zero elsewhere
VA = [c3; kron(1./(1:p)'.^2, c1./s2)];
B0 = zeros(1 + n*p, n);
B0(sub2ind(size(B0), 1 + lev, lev)) = 1;
end

function lml = soe_lml(Y, Z, nf, lev, c1, c3, s2)
% the log marginal likelihood of the VAR with block exogeneity: the sum over the two blocks
n = size(Y, 2); p = (size(Z, 2) - 1)/n; nd = n - nf;
fo = 1:nf; dm = nf+1:n;
iff = [1, reshape(1 + fo' + (0:p-1)*n, 1, [])];
[VA, B0] = soe_prior(p, n, lev, c1, c3, s2);
lml = bvar.ml.ncp(Y(:, fo), Z(:, iff), VA(iff), nf + 2, diag(s2(fo)), B0(iff, fo)) ...
    + bvar.ml.ncp(Y(:, dm), [Z, Y(:, fo)], [VA; 1./s2(fo)], n + 2, diag(s2(dm)), ...
    [B0(:, dm); zeros(nf, nd)]);
end
