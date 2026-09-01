% This function computes the graident the natural conjugate prior with
% respect to kappas
%
% See:
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019.

function [A0,VA0,nu0,S0] = prior_NCP_AD(p,kappa,Y0,Yt, three)
[Tt,n] = size(Yt);
k = 1+ n*p;
A0.v= zeros(k,n);
VA0.v = zeros(k,1);
sig2 = zeros(n,1);
if three==true
     A0.d=zeros(k*n,3);
    VA0.d=zeros(k,3);
else
    A0.d=zeros(k*n,5);
    VA0.d=zeros(k,5);
end    
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
        VA0.v(1) = kappa(3);
        VA0.d(i,3) = 1;
    else
        VA0.v(i) = kappa(1)/(l^kappa(2)*sig2(idx));
        VA0.d(i,1) = 1/(l^kappa(2)*sig2(idx));
        VA0.d(i,2) = (-log(l))*kappa(1)/(l^kappa(2)*sig2(idx));
    end
end
if three==false
    S0.v=sig2*kappa(5);
    S0.d=[zeros(n,4),sig2];
    S0=Diag(S0,n);
else
   S0=diag(sig2)*kappa(5); 
end    
nu0 = kappa(4)+n+1;
VA0 = Diag(VA0,k);
end