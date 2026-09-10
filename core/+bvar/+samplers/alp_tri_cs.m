% bvar.samplers.alp_tri_cs - Gaussian draw of the free elements of the
% unit-lower-triangular impact matrix A in the Cholesky / triangularized (CS)
% SVAR-SV model, row by row: for equation ii = 2:n the ii-1 free elements
% alp_ii regress the residual E(:,ii) on -E(:,1:ii-1) weighted by
% exp(-h(:,ii)), under independent N(0, Valp) priors. The prior mean is zero;
% there is no alp0 term.
%
%   alp = bvar.samplers.alp_tri_cs(E, h, Valp)
%   alp = bvar.samplers.alp_tri_cs(E, h, Valp, o)   % outlier-scaled (VAR-SVO)
%
%   E    : T x n residual matrix Y - XB; the CALLER computes it and the caller
%          writes the draw back into A with `A(A_id) = alp`
%   h    : T x n log-volatilities
%   Valp : n*(n-1)/2 x 1 stacked prior variances, rows ordered (2,1),
%          (3,1),(3,2), (4,1),... - row-major lower triangle
%   o    : T x 1 outlier scales, optional; default ones(T,1), which leaves the
%          weights bit-for-bit unchanged (division by 1)
%   alp  : 1 x n*(n-1)/2 ROW vector of draws; callers that want a column (as
%          ml_varsv does) transpose it themselves
%
% rng consumption: randn(ii-1,1) per equation, equations in order ii = 2:n.
%
% Provenance and the legacy copies this stands in for: tests/variant_map.md.
% Equivalence: tests/unit/test_oisv_equivalence.m,
% tests/unit/test_mlvarsv_equivalence.m.
%
% See:
% Chan, J.C.C., Koop, G. and Yu, X. (2024). Large Order-Invariant Bayesian
% VARs with Stochastic Volatility, Journal of Business and Economic
% Statistics, 42(2): 825-837.
% Chan, J.C.C. (2023). Comparing stochastic volatility specifications for large
% Bayesian VARs, Journal of Econometrics, 235(2), 1419-1446.

function alp = alp_tri_cs(E,h,Valp,o)
[T,n] = size(E);
if nargin < 4 || isempty(o)
    o = ones(T,1);
end
alp = zeros(1,n*(n-1)/2);
count_alp = 0;
for ii=2:n
    X_alpi = -E(:,1:ii-1);
    iD = sparse(1:T,1:T,exp(-h(:,ii))./o.^2);
    iValpi = sparse(1:ii-1,1:ii-1,1./Valp(count_alp+1:count_alp+ii-1));
    Kalpi = iValpi + X_alpi'*iD*X_alpi;
    alpi_hat = Kalpi\(X_alpi'*iD*E(:,ii));
    alpi = alpi_hat + chol(Kalpi,'lower')'\randn(ii-1,1);
    alp(count_alp+1:count_alp+ii-1) = alpi;
    count_alp = count_alp + ii-1;
end
end
