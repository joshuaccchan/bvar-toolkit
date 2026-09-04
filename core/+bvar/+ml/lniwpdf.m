% bvar.ml.lniwpdf - log density of the normal-inverse-Wishart distribution at
% (A,Sig): Sig ~ IW(nu0,S0), vec(A)|Sig ~ N(vec(A0), Sig kron inv(iVA0)) - the
% Kronecker-structured conjugate form, with iVA0 the PRECISION (inverse prior
% covariance) of each column's coefficients.
%
% Extracted 2026-09-02 (step 8, Kronecker family pass). Canonical source:
%   chan2020_jbes_kronecker/legacy/lniwpdf.m  (body verbatim; the only copy).
% Used by every marginal-likelihood ordinate in that package (prior ordinate
% with iVA0 = sparse diag(1./VA0), posterior ordinate with iVA0 = KA).
% Verified by tests/unit/test_kron_ml_densities.m (bitwise against the legacy
% copy) and end-to-end through the ML equivalence tests.
%
% Inputs:  A    - k x n evaluation point (VAR coefficients)
%          Sig  - n x n evaluation point (error covariance)
%          A0   - k x n mean matrix
%          iVA0 - k x k precision matrix (sparse or dense)
%          nu0  - IW degrees of freedom
%          S0   - n x n IW scale matrix
% Output:  lden - log density
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics, 38(1), 68-79.

function lden = lniwpdf(A,Sig,A0,iVA0,nu0,S0)
[k,n] = size(A);
CSig = chol(Sig,'lower');
cA = -n*k/2*log(2*pi) + n*sum(log(diag(chol(iVA0))));
cSig = -nu0*n/2*log(2) -n*(n-1)/4*log(pi) -sum(gammaln((nu0+1-(1:n))/2))...
    + nu0*sum(log(diag(chol(S0))));
tmp = A-A0;
lden = cA + cSig - (n+nu0+k+1)*sum(log(diag(CSig))) ...
    - .5*trace(Sig\(S0+tmp'*iVA0*tmp));
end
