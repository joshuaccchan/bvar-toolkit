% bvar.samplers.factor_fsv - joint Gaussian draw of the whole latent factor path
% F (T x r) in the factor-SV VAR, conditional on the loadings L, the VAR
% coefficients A and the n+r log-volatilities h. Precision sampler on the
% stacked T*r vector; the idiosyncratic variances are exp(h(:,1:n)), the factor
% variances exp(h(:,n+1:n+r)).
%
%   F = bvar.samplers.factor_fsv(Y, X, A, L, h)
%
%   L : n x r loading matrix (unit-lower-triangular top block)
%   h : T x (n+r), idiosyncratic columns first
%
% See:
% Chan, J.C.C. (2023). Comparing stochastic volatility specifications for large
% Bayesian VARs, Journal of Econometrics, 235(2), 1419-1446.

function F = factor_fsv(Y,X,A,L,h)
[T,n] = size(Y);
r = size(L,2);
e = reshape((Y-X*A)',T*n,1);
Xf = kron(speye(T),L);
XfiSig = Xf'*sparse(1:T*n,1:T*n,reshape(exp(-h(:,1:n))',T*n,1));
Kf = sparse(1:T*r,1:T*r,reshape(exp(-h(:,n+1:end))',T*r,1)) + XfiSig*Xf;
f_hat = Kf\(XfiSig*e);
f = f_hat + chol(Kf,'lower')'\randn(T*r,1);
F = reshape(f,r,T)';
end
