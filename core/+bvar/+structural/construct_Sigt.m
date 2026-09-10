% bvar.structural.construct_Sigt - time-varying reduced-form covariance matrices
% Sigt(t,:,:) = B0^{-1} diag(exp(h_t)) B0^{-T} implied by an impact matrix B0
% (full n x n in the OI model; unit-lower-triangular A in the CS model) and the
% T x n log-volatility paths h. Pure transform - consumes no rng.
%
% Provenance and the legacy copies this stands in for: tests/variant_map.md.
% Equivalence: tests/unit/test_construct_sigt.m.
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
