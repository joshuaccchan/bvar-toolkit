% bvar.samplers.eq_tri_cs - equation-by-equation Gaussian draw of the VAR
% coefficient matrix B (n x k rows, intercept first) in the CHOLESKY /
% triangularized (CS) SVAR-SV model, conditioning on the unit-lower-triangular
% impact matrix A and the log-volatilities h. Equation ii stacks the rows
% ii:n of the rotated system (Ytilde = Y*A', partialling out the other
% equations' fits through XB*A') and draws B(ii,:) from its Gaussian
% conditional; XB(:,ii) = X*betai is refreshed in place so later equations
% condition on the new draw. Prior variances enter as the stacked k*n vector
% tmpdV, which the caller supplies (the Vbeta output of bvar.priors.vtheta
% produces it); prior means as the stacked k*n vector beta0.
%
%   [B,XB] = bvar.samplers.eq_tri_cs(Y, X, XB, B, A, h, tmpdV, beta0)
%
% THIS IS THE CORRECTED TRIANGULAR ALGORITHM of Carriero, Chan, Clark and
% Marcellino (2022), the corrigendum to Carriero, Clark and Marcellino (2019).
% The original algorithm drew equation j from a conditional that omitted part of
% the information - it conditioned on y(1),...,y(j-1) rather than on the whole of
% y - so it did not sample the intended triangular factorization. The corrigendum
% keeps that factorization and restores the missing term at the same O(n^4) cost.
% Stacking rows ii:n above, rather than equation ii alone, is that correction.
%
% See:
% Carriero, A., Chan, J.C.C., Clark, T.E. and Marcellino, M. (2022). Corrigendum
% to: Large Bayesian Vector Autoregressions with Stochastic Volatility and
% Non-Conjugate Priors, Journal of Econometrics, 227(2): 506-512.
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
