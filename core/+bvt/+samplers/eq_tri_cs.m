% bvt.samplers.eq_tri_cs - equation-by-equation Gaussian draw of the VAR
% coefficient matrix B (n x k rows, intercept first) in the CHOLESKY /
% triangularized (CS) SVAR-SV model, conditioning on the unit-lower-triangular
% impact matrix A and the log-volatilities h. Equation ii stacks the rows
% ii:n of the rotated system (Ytilde = Y*A', partialling out the other
% equations' fits through XB*A') and draws B(ii,:) from its Gaussian
% conditional; XB(:,ii) = X*betai is refreshed in place so later equations
% condition on the new draw. Prior variances enter as the stacked k*n vector
% tmpdV (caller-computed - legacy getVbeta(idx_kappa1,idx_kappa2,kappa,
% C.*Psi,sig2), reproduced exactly by the Vbeta output of bvt.priors.vtheta);
% prior means as the stacked k*n vector beta0 (legacy Hyper.beta0).
%
%   [B,XB] = bvt.samplers.eq_tri_cs(Y, X, XB, B, A, h, tmpdV, beta0)
%
% rng consumption: randn(k,1) per equation, equations in order ii = 1:n.
% The caller keeps `beta = reshape(B',k_beta,1)` (legacy line 73).
%
% Extracted 2026-09-02 (step 7, OISV family pass). Canonical source (body
% verbatim): chan_koop_yu2024_jbes_oisv/legacy/CS_MH.m lines 54-72 (the inline
% "sample B" block, from the Ytilde line through the equation loop). Also
% canonicalizes forecast_CS_MH.m lines 45-63, textually identical modulo the
% Y/X/T -> Yt/Xt/Tt renaming (all enter through the arguments/sizes here).
% Edits made, in full: wrapped as a function with [T,n] = size(Y) and
% k = size(X,2) replacing the workspace T, n, k = k_beta/n (identical
% integers); Hyper.beta0 renamed beta0. Everything else byte-verbatim.
% Draw-for-draw equivalence: tests/unit/test_oisv_equivalence.m.
%
% See:
% Chan, J.C.C., Koop, G. and Yu, X. (2024). Large Order-Invariant Bayesian
% VARs with Stochastic Volatility, Journal of Business and Economic
% Statistics, 42(2): 825-837.

function [B,XB] = eq_tri_cs(Y,X,XB,B,A,h,tmpdV,beta0)
[T,n] = size(Y);
k = size(X,2);
Ytilde = Y*sparse(A');
for ii = 1:n
    tmpXBA = XB(:,[1:ii-1 ii+1:n])*sparse(A(:,[1:ii-1 ii+1:n])');
    Zi = Ytilde(:,ii:n) - tmpXBA(:,ii:n);
    zi = reshape(Zi',T*(n-ii+1),1); %#ok<NASGU> % kept verbatim from the legacy block (assigned there and unused there too)
    Xi = repmat(X,n-ii+1,1);
    tmp1 = exp(-h(:,ii:n)).*repmat(A(ii:n,ii)',T,1);
    tmp2 = tmp1.*repmat(A(ii:n,ii)',T,1);
    iVbetai = sparse(1:k,1:k,1./tmpdV((ii-1)*k+1:ii*k));
    betai0 = beta0((ii-1)*k+1:ii*k);
    Kbetai = iVbetai + Xi'*sparse(1:(n-ii+1)*T,1:(n-ii+1)*T,tmp2(:))*Xi;
    CKbetai = chol(Kbetai,'lower');
    betai_hat = (CKbetai')\(CKbetai\(iVbetai*betai0 ...
        + Xi'*sparse(1:(n-ii+1)*T,1:(n-ii+1)*T,tmp1(:))*Zi(:)));

    betai = betai_hat + CKbetai'\randn(k,1);
    B(ii,:) = betai;
    XB(:,ii) = X*betai;
end
end
