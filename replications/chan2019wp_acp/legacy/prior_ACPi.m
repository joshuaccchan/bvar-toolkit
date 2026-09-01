% This script constructs the asymmetric conjugate prior given the
% hyperparameter values
%
% See:
% Chan, J.C.C. (2019). Asymmetric conjugate priors for large Bayesian VARs,
% CAMA Working Papers 51/2019

function [mi,Vi,nui,Si] = prior_ACPi(n,p,var_i,kappa,sig2)
ki = var_i + n*p;
mi = zeros(ki,1);
Vi = zeros(ki,1);
    % construct Vi
for j=1:ki
    if j <= n*p+1
        l = ceil((j-1)/n); % lag length
        idx = mod(j-1,n);  % variable index
        if idx==0
            idx = n;
        end    
    else    
        idx = j - (n*p+1);
    end
    
    if j==1 % intercept
        Vi(j) = kappa(4);
    elseif j > n*p+1    % alpha_i
        Vi(j) = kappa(3)/sig2(idx);
    elseif idx == var_i % own lag
        Vi(j) = kappa(1)/(l^2*sig2(idx));
    else % lag of other variables
        Vi(j) = kappa(2)/(l^2*sig2(idx));
    end    
end
Si = sig2(var_i)/2; 
nui = 1 + var_i/2;
Vi = sparse(1:ki,1:ki,Vi);
end