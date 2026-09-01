% This function constructs the natural conjugate prior
%
% See:
% Chan, J.C.C. (2022). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, forthcoming.

function [A0,VA0,nu0,S0,U_hat] = prior_NCP(p,c1,c2,Y0,Y)
[T,n] = size(Y);
k = 1+n*p;
A0 = zeros(k,n);
VA0 = zeros(k,1);
sig2 = zeros(n,1);
U_hat = zeros(T,n);
    % construct VA0
tmpY = [Y0(end-4+1:end,:); Y];
for i=1:n
    Z = [ones(T,1) tmpY(4:end-1,i) tmpY(3:end-2,i) tmpY(2:end-3,i)...
        tmpY(1:end-4,i)];
    tmpb = (Z'*Z)\(Z'*tmpY(5:end,i));
    U_hat(:,i) = tmpY(5:end,i)-Z*tmpb;
    sig2(i) = mean(U_hat(:,i).^2);
end
for i=1:k
    l = ceil((i-1)/n);
    idx = mod(i-1,n); % variable index
    if idx==0
        idx = n;
    end
    if i==1 % intercept
        VA0(1) = c2;
    else
        VA0(i) = c1/(l^2*sig2(idx));
    end
end
S0 = diag(sig2); nu0 = n+3;
end