% bvar.util.surform2 - sparse SUR expansion with rows kron(speye(n), X(t,:)) for
% stacked-vector VAR sampling.
%
%   Xout = bvar.util.surform2(X, n)
%
%   X    : T x k matrix of regressors
%   n    : number of equations
%   Xout : T*n x n*k sparse matrix whose tth block of n rows is
%          kron(speye(n), X(t,:)), so that Xout*vec(A) returns the fitted
%          values stacked as reshape(Y',T*n,1), where A is the k x n
%          coefficient matrix
%
% NOT the same operator as bvar.util.surform, which is the T x Tk
% block-diagonal expansion used for TVP state stacking.
%
% Provenance and the legacy copies this stands in for: tests/variant_map.md.
% Equivalence: tests/unit/test_surform2.m.
%
% See:
% Chan, J.C.C. (2023). Comparing stochastic volatility specifications for large
% Bayesian VARs, Journal of Econometrics, 235(2), 1419-1446.
function Xout = surform2( X, n )
repX = kron(X,ones(n,1));
[r,c] = size( X );
idi = kron((1:r*n)',ones(c,1));
idj = repmat((1:n*c)',r,1);
Xout = sparse(idi,idj,reshape(repX',n*r*c,1));
end