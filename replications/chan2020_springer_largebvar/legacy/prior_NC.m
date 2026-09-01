% This script constructs the natural conjugate prior
%
% See:
% Chan, J.C.C. (2020). Large Bayesian Vector Autoregressions. In: P. Fuleky (Eds),
% Macroeconomic Forecasting in the Era of Big Data, 95-125, Springer, Cham

function [A0,VA0,nu0,S0] = prior_NC(p,c1,c2,Y0,shortYt)
[Tt,n] = size(shortYt);
k = 1+ n*p;
A0 = zeros(k,n);
VA0 = zeros(k,1);
sig2 = zeros(n,1);
    % construct VA0
tmpY = [Y0(end-p+1:end,:); shortYt];
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
        VA0(1) = c2;
    else
        VA0(i) = c1/(l^2*sig2(idx));
    end
end
S0 = diag(sig2); nu0 = n+3;
end