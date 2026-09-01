% Extracted 2026-09-01 (step 3, zero-risk core): verified identical (modulo comments/whitespace). Canonical source: chan2023_joe_mlvarsv/legacy/utility/ldet.m (also code-identical in cjz2019_ad_opthyper/legacy).
% This function evaluates the log determinant
%
% See:
% Chan, J.C.C. (2022). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, forthcoming.

function k = ldet(Omega)
    k = 2*sum(log(diag(chol(Omega,'lower'))));
end