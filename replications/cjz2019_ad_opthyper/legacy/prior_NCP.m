% This function constructs the natural conjugate prior
%
% See:
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019.

function [A0,VA0,nu0,S0] = prior_NCP(p,kappa,Y0,Yt)
[Tt,n] = size(Yt);
k = 1+ n*p;
A0 = zeros(k,n);
VA0 = zeros(k,1);
sig2 = zeros(n,1);
    % construct VA0
tmpY = [Y0(end-p+1:end,:); Yt];
for i=1:n
    Z = [ones(Tt,1) tmpY(4:end-1,i) tmpY(3:end-2,i) tmpY(2:end-3,i)...
        tmpY(1:end-4,i)];
    tmpb = (Z'*Z)\(Z'*tmpY(5:end,i));
    sig2(i) = mean((tmpY(5:end,i)-Z*tmpb).^2);
end
for i=1:k
    l = ceil((i-1)/n);
    idx = mod(i-1,n); % variable index
    if idx==0
        idx = n;
    end
    if i==1 % intercept
        VA0(1) = kappa(3);
    else
        VA0(i) = kappa(1)/(l^kappa(2)*sig2(idx));
    end
end
S0 = kappa(5)*diag(sig2); nu0 = kappa(4)+n+1;
VA0 = sparse(1:k,1:k,VA0);
end