% Body from chan2023_joe_mlvarsv/legacy/utility/ldet.m (identical modulo comments).
% Equivalence: tests/unit/test_util_small.m. Record: tests/variant_map.md.
% This function evaluates the log determinant
%
% See:
% Chan, J.C.C. (2023). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, 235(2): 1419-1446.

function k = ldet(Omega)
    k = 2*sum(log(diag(chol(Omega,'lower'))));
end