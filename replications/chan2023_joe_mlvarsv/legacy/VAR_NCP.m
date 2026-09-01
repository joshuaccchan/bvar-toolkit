% This script estimates the standard VAR with the natural conjugate prior
% and computes the associated marginal likelihood
% (MCMC is not needed)
%
% See:
% Chan, J.C.C. (2022). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, forthcoming.

% compute the parameters for the posterior density
XX = X'*X;
A_tilde = XX\(X'*Y);
K_A = sparse(1:k,1:k,1./Hyper.VA) + XX;
    % posterior mean of the VAR coefficients, arranged as a k by n matrix
A_hat = K_A\(sparse(1:k,1:k,Hyper.VA)\Hyper.A0 + XX*A_tilde); 
S_hat = Hyper.S0 + Hyper.A0'*sparse(1:k,1:k,1./Hyper.VA)*Hyper.A0 + Y'*Y - A_hat'*K_A*A_hat;
S_hat = (S_hat+S_hat')/2;
    
if cp_ml    
    lml = -n*T/2*log(pi) - n/2*(sum(log(Hyper.VA)) + ldet(K_A)) + Hyper.nu0/2*ldet(Hyper.S0)...
        - (Hyper.nu0+T)/2*ldet(S_hat) + mgammaln(n,(Hyper.nu0+T)/2) - mgammaln(n,Hyper.nu0/2);    
end
