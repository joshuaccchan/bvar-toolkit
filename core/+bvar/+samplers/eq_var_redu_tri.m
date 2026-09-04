% bvar.samplers.eq_var_redu_tri - equation-by-equation Gaussian draw of the VAR
% coefficient matrix A (k x n, intercept first) in the reduced-form VAR-SV with
% a unit-lower-triangular impact matrix B0. Equation ii zeroes A(:,ii), stacks
% only rows ii:n of the rotated system (rows 1..ii-1 have B0(j,ii) = 0, so they
% contribute nothing), and draws A(:,ii) from its Gaussian conditional.
%
%   A = bvar.samplers.eq_var_redu_tri(Y, X, B0, h, A, Valp, alp0)      % VAR-SV
%   A = bvar.samplers.eq_var_redu_tri(Y, X, B0, h, A, Valp, alp0, o)   % VAR-SVO
%
%   Valp : k*n x 1 stacked prior variances (legacy Hyper.Valp)
%   alp0 : k*n x 1 stacked prior means (legacy Hyper.alp0)
%   o    : T x 1 outlier scales, optional; default ones(T,1), which leaves the
%          row scaling bit-for-bit unchanged (multiplication by 1)
%
% rng: randn(k,1) per equation, ii = 1:n. The caller keeps alp = reshape(A,k_alp,1).
%
% Sibling of bvar.samplers.eq_svar_oi (order-invariant SVAR-SV), which stacks all
% n rows, has no prior-mean term and no o. Do not merge the two.
%
% Extracted 2026-09-03 (step 9). Canonical source, body verbatim:
% chan2023_joe_mlvarsv/legacy/VAR_ARSV_redu.m lines 44-57; the o argument also
% canonicalizes VAR_ARSVO_redu.m lines 51-64 (sole textual difference: the
% Lambda line's .*repmat(o,1,n-ii+1)). Wrapped as a function with sizes from the
% arguments and vec -> bvar.util.vec. Details: tests/variant_map.md.
%
% See:
% Chan, J.C.C. (2023). Comparing stochastic volatility specifications for large
% Bayesian VARs, Journal of Econometrics, 235(2), 1419-1446.

function A = eq_var_redu_tri(Y,X,B0,h,A,Valp,alp0,o)
[T,n] = size(Y);
k = size(X,2);
if nargin < 8 || isempty(o)
    o = ones(T,1);
end
sqrt_exph = exp(h/2);
for ii = 1:n
    A(:,ii) = 0;
    Lambda = bvar.util.vec(sqrt_exph(:,ii:n).*repmat(o,1,n-ii+1));
    yi = bvar.util.vec((Y-X*A)*B0(ii:n,:)')./Lambda;
    Wi = kron(B0(ii:n,ii),X)./Lambda;
    iValpi = sparse(1:k,1:k,1./Valp((ii-1)*k+1:ii*k));
    alpi0 = alp0((ii-1)*k+1:ii*k);
    Kalpi = iValpi + Wi'*Wi;
    CKalpi = chol(Kalpi,'lower');
    alpi_hat = (CKalpi')\(CKalpi\(iValpi*alpi0 + Wi'*yi));
    alpi = alpi_hat + CKalpi'\randn(k,1);
    A(:,ii) = alpi;
end
end
