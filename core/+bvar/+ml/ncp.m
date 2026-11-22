% bvar.ml.ncp - log marginal likelihood of the multivariate regression Y = X*B + E under the
% natural conjugate prior, in closed form, and the posterior. Written for this toolkit.
%
%   [lml, Bhat, K, Shat] = bvar.ml.ncp(Y, X, VA, nu0, S0)
%   [lml, Bhat, K, Shat] = bvar.ml.ncp(Y, X, VA, nu0, S0, B0)
%
%   Y    : T x n observations; a VAR has the rows y_t'
%   X    : T x k regressors; a VAR has the rows (1, y_{t-1}', ..., y_{t-p}') of bvar.util.build_lags
%   VA   : k x 1 prior variances, positive: vec(B) | Sig ~ N(vec(B0), Sig kron diag(VA))
%   nu0  : prior degrees of freedom, Sig ~ IW(nu0, S0), with nu0 > n - 1
%   S0   : n x n prior scale matrix, positive definite
%   B0   : k x n prior mean, default zeros
%   lml  : the log marginal likelihood
%   Bhat, K, Shat : the posterior, vec(B) | Sig, Y ~ N(vec(Bhat), Sig kron inv(K)) and
%          Sig | Y ~ IW(nu0 + T, Shat)
%
% See:
% Chan, J.C.C. (2027). Bayesian Macroeconometrics: Methods and Applications. Chapman &
% Hall/CRC.

function [lml, Bhat, K, Shat] = ncp(Y, X, VA, nu0, S0, B0)
[T, n] = size(Y);
k = size(X, 2);
if size(X, 1) ~= T
    error('bvar:ml:ncp:badDims', 'X has %d rows and Y has %d', size(X, 1), T);
end
if numel(VA) ~= k || any(VA(:) <= 0)
    error('bvar:ml:ncp:badPrior', 'VA must contain %d positive prior variances', k);
end
if ~isequal(size(S0), [n n]) || ~(isscalar(nu0) && nu0 > n - 1)
    error('bvar:ml:ncp:badPrior', 'S0 must be %d x %d and nu0 greater than %d', n, n, n - 1);
end
if nargin < 6
    B0 = zeros(k, n);
elseif ~isequal(size(B0), [k n])
    error('bvar:ml:ncp:badPrior', 'B0 must be %d x %d', k, n);
end
VA = VA(:);
K = diag(1./VA) + X'*X;
Bhat = K\(B0./VA + X'*Y);
E = Y - X*Bhat; D = Bhat - B0;
Shat = S0 + E'*E + D'*(D./VA);           % S0 + Y'*Y + B0'*(B0./VA) - Bhat'*K*Bhat, computed stably
Shat = (Shat + Shat')/2;
lml = -n*T/2*log(pi) - n/2*(sum(log(VA)) + bvar.util.ldet(K)) ...
    + nu0/2*bvar.util.ldet(S0) - (nu0 + T)/2*bvar.util.ldet(Shat) ...
    + bvar.util.mgammaln(n, (nu0 + T)/2) - bvar.util.mgammaln(n, nu0/2);
end
