% bvt.structural.construct_Sigt - time-varying reduced-form covariance matrices
% Sigt(t,:,:) = B0^{-1} diag(exp(h_t)) B0^{-T} implied by an impact matrix B0
% (full n x n in the OI model; unit-lower-triangular A in the CS model) and the
% T x n log-volatility paths h. Pure transform - consumes no rng.
%
% Extracted 2026-09-02 (step 7, OISV family pass). Canonical source (body
% verbatim): chan_koop_yu2024_jbes_oisv/legacy/utility/construct_Sigt.m.
% Also canonicalizes the private subfunction copy inside
% chan_koop_yu2024_jbes_oisv/legacy/func_main_SVAR_v2.m (lines 67-73), which is
% comment-stripped IDENTICAL (verified by diff 2026-09-02; at runtime the legacy
% func resolves its private copy, the harness/utility path resolves the utility
% copy - same code either way). Only the namespace was added; nothing renamed
% or parameterized.
%
% This function constructs the time-varying covariance matrices
%
% See:
% Chan, J.C.C., Koop, G. and Yu, X. (2024). Large Order-Invariant Bayesian
% VARs with Stochastic Volatility, Journal of Business and Economic
% Statistics, 42(2): 825-837.

function Sigt = construct_Sigt(h,B0)
[T,n] = size(h);
Sigt = zeros(T,n,n);
for t=1:T
    Sigt(t,:,:) = (B0\sparse(1:n,1:n,exp(h(t,:))))/(B0');
end
end
