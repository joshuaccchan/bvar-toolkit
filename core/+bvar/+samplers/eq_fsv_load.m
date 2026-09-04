% bvar.samplers.eq_fsv_load - equation-by-equation Gaussian draw of the VAR
% coefficients A (k x n) jointly with the free factor loadings L (n x r) in the
% factor-SV VAR, conditional on the factor path F and the log-volatilities h.
% Equation ii regresses Y(:,ii) on [X F(:,1:min(ii-1,r))] weighted by
% exp(-h(:,ii)); the first r equations are normalized (unit diagonal, so
% F(:,ii) is subtracted from the left-hand side instead of loaded).
%
%   [A,L] = bvar.samplers.eq_fsv_load(Y, X, F, h, A, L, Valp, alp0, Vl, l0)
%
%   Valp, alp0 : k*n x 1 stacked prior variances / means of A
%   Vl, l0     : scalar prior variance / mean of every free loading
%
% rng: randn(k+min(ii-1,r),1) per equation, ii = 1:n. The caller keeps alp = A(:).
%
% Extracted 2026-09-03 (step 9). Canonical source, body verbatim:
% chan2023_joe_mlvarsv/legacy/VAR_FSV.m lines 48-71 (single copy in the repo).
% Wrapped as a function with T, n, k, r taken from the arguments and Hyper.*
% passed explicitly.
%
% See:
% Chan, J.C.C. (2023). Comparing stochastic volatility specifications for large
% Bayesian VARs, Journal of Econometrics, 235(2), 1419-1446.

function [A,L] = eq_fsv_load(Y,X,F,h,A,L,Valp,alp0,Vl,l0)
[T,n] = size(Y);
k = size(X,2);
r = size(F,2);
for ii = 1:n
    if ii <= r
        Zi = [X F(:,1:ii-1)];
        iVthetai = sparse(1:k+ii-1,1:k+ii-1,1./[Valp((ii-1)*k+1:ii*k); ...
            Vl*ones(ii-1,1)]);
        thetai0 = [alp0((ii-1)*k+1:ii*k); l0*ones(ii-1,1)];
        ZiSigi = Zi'*sparse(1:T,1:T,exp(-h(:,ii)));
        dthetai = iVthetai*thetai0 + ZiSigi*(Y(:,ii)-F(:,ii));
    else
        Zi = [X F];
        iVthetai = sparse(1:k+r,1:k+r,1./[Valp((ii-1)*k+1:ii*k); ...
            Vl*ones(r,1)]);
        thetai0 = [alp0((ii-1)*k+1:ii*k); l0*ones(r,1)];
        ZiSigi = Zi'*sparse(1:T,1:T,exp(-h(:,ii)));
        dthetai = iVthetai*thetai0 + ZiSigi*Y(:,ii);
    end
    Kthetai = iVthetai + ZiSigi*Zi;
    CKthetai = chol(Kthetai,'lower');
    thetai_hat = CKthetai'\(CKthetai\dthetai);
    thetai = thetai_hat + CKthetai'\randn(k+min(ii-1,r),1);

    A(:,ii) = thetai(1:k);
    L(ii,1:min(ii-1,r)) = thetai(k+1:end);
end
end
