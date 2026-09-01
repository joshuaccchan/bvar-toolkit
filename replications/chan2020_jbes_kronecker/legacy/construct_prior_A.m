% This script constructs the prior for the VAR coefficients
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

A0 = zeros(k,n);
VA0 = zeros(k,1);
sig2 = zeros(n,1);
tmpY = [Y0(end-p+1:end,:); shortY];
c1 = .2^2; c2 = 100; % hyperparameters for slopes and intercepts
for i=1:n
    Z = [ones(T,1) tmpY(4:end-1,i) tmpY(3:end-2,i) tmpY(2:end-3,i)...
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