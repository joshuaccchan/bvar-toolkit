% bvt.samplers.alp_tri_cs - Gaussian draw of the free elements of the
% unit-lower-triangular impact matrix A in the Cholesky / triangularized (CS)
% SVAR-SV model, row by row: for equation ii = 2:n the ii-1 free elements
% alp_ii regress the residual E(:,ii) on -E(:,1:ii-1) weighted by
% exp(-h(:,ii)), with independent N(0, Valp) priors (zero prior mean - the
% legacy block never adds an alp0 term).
%
%   alp = bvt.samplers.alp_tri_cs(E, h, Valp)
%
%   E    : T x n residual matrix Y - XB (computed by the CALLER, verbatim
%          legacy position; the caller also keeps `A(A_id) = alp`)
%   h    : T x n log-volatilities
%   Valp : n*(n-1)/2 x 1 stacked prior variances (legacy Hyper.Valp), rows
%          ordered (2,1), (3,1),(3,2), (4,1),... - row-major lower triangle
%   alp  : 1 x n*(n-1)/2 row vector of draws (the legacy scripts grow `alp`
%          dynamically into exactly this 1 x k_alp row in the first sweep and
%          fully overwrite it every sweep; preallocated fresh here,
%          value-identical)
%
% rng consumption: randn(ii-1,1) per equation, equations in order ii = 2:n.
%
% Extracted 2026-09-02 (step 7, OISV family pass). Canonical source (body
% verbatim): chan_koop_yu2024_jbes_oisv/legacy/CS_MH.m lines 77-87 (the inline
% "sample alp" count_alp loop). Also canonicalizes forecast_CS_MH.m lines
% 68-78, textually identical modulo T -> Tt (enters through size(E) here).
% Edits made, in full: wrapped as a function with [T,n] = size(E) replacing
% the workspace T, n (identical integers); Hyper.Valp renamed Valp; alp
% preallocated zeros(1,n*(n-1)/2) instead of dynamically grown (see above).
% Everything else byte-verbatim. Draw-for-draw equivalence:
% tests/unit/test_oisv_equivalence.m.
%
% See:
% Chan, J.C.C., Koop, G. and Yu, X. (2024). Large Order-Invariant Bayesian
% VARs with Stochastic Volatility, Journal of Business and Economic
% Statistics, 42(2): 825-837.

function alp = alp_tri_cs(E,h,Valp)
[T,n] = size(E);
alp = zeros(1,n*(n-1)/2);
count_alp = 0;
for ii=2:n
    X_alpi = -E(:,1:ii-1);
    iD = sparse(1:T,1:T,exp(-h(:,ii)));
    iValpi = sparse(1:ii-1,1:ii-1,1./Valp(count_alp+1:count_alp+ii-1));
    Kalpi = iValpi + X_alpi'*iD*X_alpi;
    alpi_hat = Kalpi\(X_alpi'*iD*E(:,ii));
    alpi = alpi_hat + chol(Kalpi,'lower')'\randn(ii-1,1);
    alp(count_alp+1:count_alp+ii-1) = alpi;
    count_alp = count_alp + ii-1;
end
end
