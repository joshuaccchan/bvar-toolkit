% bvar.samplers.eq_svar_oi - equation-by-equation Gaussian draw of the VAR
% coefficient matrix A (k x n, intercept first) in the ORDER-INVARIANT
% structural SVAR-SV model: for each equation ii the column A(:,ii) is zeroed,
% the standardized system yi = vec((Y-X*A)*B0')./Lambda with regressors
% kron(B0(:,ii),X)./Lambda (Lambda = vec(exp(h/2))) is formed, and
% A(:,ii) is drawn from its Gaussian conditional. Prior variances enter as the
% stacked k*n vector tmpdV (the caller computes it - legacy
% getVbeta(idx_kappa1,idx_kappa2,kappa,C.*Psi,sig2), reproduced exactly by the
% Vbeta output of bvar.priors.vtheta on the same inputs).
%
%   A = bvar.samplers.eq_svar_oi(Y, X, B0, h, A, tmpdV)
%
% rng consumption: randn(k,1) per equation, equations in order ii = 1:n.
% The caller keeps `alpha = A(:)` (legacy line 88).
%
% Body from chan_koop_yu2024_jbes_oisv/legacy/SVARSV_MH.m lines 76-87 (the
% inline "sample alpha" block), wrapped as a function: k and n from the
% arguments, vec -> bvar.util.vec. It does NOT cover forecast_SVARSV_MH.m lines
% 69-81, which rewrites this step with time-interleaved stacking
% (Wi = kron(Xt,B0(:,ii))) and an explicit exp(-h) weighting matrix in place of
% the ./Lambda row scaling. That block draws from the same conditional but orders
% the floating-point operations differently, so its draws are not bitwise equal
% to these.
% Equivalence: tests/unit/test_oisv_equivalence.m. Record: tests/variant_map.md.
%
% See:
% Chan, J.C.C., Koop, G. and Yu, X. (2024). Large Order-Invariant Bayesian
% VARs with Stochastic Volatility, Journal of Business and Economic
% Statistics, 42(2): 825-837.

function A = eq_svar_oi(Y,X,B0,h,A,tmpdV)
k = size(X,2);
n = size(Y,2);
Lambda = bvar.util.vec(exp(h/2));
for ii=1:n
    A(:,ii) = 0;
    yi = bvar.util.vec((Y-X*A)*B0')./Lambda;
    Wi = kron(B0(:,ii),X)./Lambda;
    iValphai = sparse(1:k,1:k,1./tmpdV((ii-1)*k+1:ii*k));
    Kalphai = iValphai + Wi'*Wi;
    CKalphai = chol(Kalphai,'lower');
    alphai_hat = (CKalphai')\(CKalphai\(Wi'*yi));
    alphai = alphai_hat + CKalphai'\randn(k,1);
    A(:,ii) = alphai;
end
end
