% bvt.ml.linvgammpdf - log density of the inverse-gamma IG(a,b) distribution
% at y (shape a, scale b; elementwise over array inputs).
%
% Extracted 2026-09-02 (step 8, Kronecker family pass). Canonical source:
%   chan2020_jbes_kronecker/legacy/linvgammpdf.m  (body verbatim; the only
%   copy). Used by the sigh2 posterior ordinates of the CSV-family marginal
% likelihoods. Verified by tests/unit/test_kron_ml_densities.m (bitwise) and
% end-to-end through the ML equivalence tests.
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics, 38(1), 68-79.

function lden = linvgammpdf(y, a, b)

lden = a.*log(b) - gammaln(a) - (a+1) .* log(y) - b./y;

end
