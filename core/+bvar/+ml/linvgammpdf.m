% bvar.ml.linvgammpdf - log density of the inverse-gamma IG(a,b) distribution
% at y (shape a, scale b; elementwise over array inputs).
%
% Body from chan2020_jbes_kronecker/legacy/linvgammpdf.m.
% Equivalence: tests/unit/test_kron_ml_densities.m. Record: tests/variant_map.md.
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics, 38(1), 68-79.

function lden = linvgammpdf(y, a, b)

lden = a.*log(b) - gammaln(a) - (a+1) .* log(y) - b./y;

end
