% bvar.samplers.alp_tri_cs - Gaussian draw of the free elements of the
% unit-lower-triangular impact matrix A in the Cholesky / triangularized (CS)
% SVAR-SV model, row by row: for equation ii = 2:n the ii-1 free elements
% alp_ii regress the residual E(:,ii) on -E(:,1:ii-1) weighted by
% exp(-h(:,ii)), with independent N(0, Valp) priors (zero prior mean - the
% legacy block never adds an alp0 term).
%
%   alp = bvar.samplers.alp_tri_cs(E, h, Valp)
%   alp = bvar.samplers.alp_tri_cs(E, h, Valp, o)   % outlier-scaled (VAR-SVO)
%
%   E    : T x n residual matrix Y - XB (computed by the CALLER, verbatim
%          legacy position; the caller also keeps `A(A_id) = alp`)
%   h    : T x n log-volatilities
%   Valp : n*(n-1)/2 x 1 stacked prior variances (legacy Hyper.Valp), rows
%          ordered (2,1), (3,1),(3,2), (4,1),... - row-major lower triangle
%   o    : T x 1 outlier scales, optional; default ones(T,1), which leaves the
%          weights bit-for-bit unchanged (division by 1)
%   alp  : 1 x n*(n-1)/2 row vector of draws (the legacy scripts grow `alp`
%          dynamically into exactly this 1 x k_alp row in the first sweep and
%          fully overwrite it every sweep; preallocated fresh here,
%          value-identical). ml_varsv keeps the same draw as a column - its
%          caller transposes.
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
% Also canonicalizes (added 2026-09-03, step 9; draw-for-draw verified, not a
% visual match): the impact-matrix block of chan2023_joe_mlvarsv, i.e.
% VAR_ARSV_redu.m lines 64-73 (o omitted) and VAR_ARSVO_redu.m lines 71-80
% (o supplied - its sole textual difference is the iD line's ./o.^2). Those
% callers name the object beta/Hyper.Vbeta and keep it as a k_beta x 1 column,
% so they transpose the returned row. Details: tests/variant_map.md.
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
