% bvar.ml.llike_ma - log likelihood of the homoskedastic BVAR-MA(1) model at
% MA coefficient psi, given the untransformed residuals U = Y - X*A and error
% covariance Sig: u_t = eps_t + psi*eps_{t-1}, eps_t ~ N(0,Sig), with the
% first transformed observation carrying the initialization variance factor
% (1+psi^2) and the matching -n/2*log(1+psi^2) constant.
%
% Extracted 2026-09-02 (step 8, Kronecker family pass). Canonical source:
%   chan2020_jbes_kronecker/legacy/llike_MA.m  (body verbatim).
% The realtime_forecasts/llike_MA.m file in the same package is NOT
% canonicalized here: its function line is named llike_MA1 (filename wins at
% dispatch) and it differs only in comments/whitespace - it belongs to the
% part-2 forecast pass. Note the legacy body writes CSig = chol(Sig)' (upper
% Cholesky transposed) where the rest of the package writes chol(Sig,'lower');
% kept verbatim. Verified by tests/unit/test_kron_ml_densities.m (bitwise)
% and end-to-end through the BVAR-MA equivalence test.
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics, 38(1), 68-79.

function llike = llike_ma(psi,U,Sig)
[T,n] = size(U);
Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
Utld = Hpsi\U;
CSig = chol(Sig)';
c = -T*n/2*log(2*pi) - T*sum(log(diag(CSig))) - n/2*log(1+psi^2);
tmp = (Utld/CSig');
s2 = sum(tmp.^2,2);
llike = c -.5*(s2(1)/(1+psi^2) + sum(s2(2:end)));
end
