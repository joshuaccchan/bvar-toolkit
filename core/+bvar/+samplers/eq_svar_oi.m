% bvar.samplers.eq_svar_oi - equation-by-equation Gaussian draw of the VAR
% coefficient matrix A (k x n, intercept first) in the ORDER-INVARIANT
% structural SVAR-SV model: for each equation ii the column A(:,ii) is zeroed,
% the standardized system yi = vec((Y-X*A)*B0')./Lambda with regressors
% kron(B0(:,ii),X)./Lambda (Lambda = vec(exp(h/2))) is formed, and
% A(:,ii) is drawn from its Gaussian conditional. Prior variances enter as the
% stacked k*n vector tmpdV, which the caller supplies; the Vbeta output of
% bvar.priors.vtheta produces it.
%
%   A = bvar.samplers.eq_svar_oi(Y, X, B0, h, A, tmpdV)
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
