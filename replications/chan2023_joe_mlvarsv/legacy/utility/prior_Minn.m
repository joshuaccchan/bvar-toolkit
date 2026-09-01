% This script constructs the Minnesota prior
%
% See:
% Chan, J.C.C. (2022). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, forthcoming.

function [beta_Minn,V_Minn,Sig_hat,U_hat] = prior_Minn(p,c1,c2,c3,Y0,Yt)
[Tt,n] = size(Yt);
k = 1 + n*p;
beta_Minn = zeros(k*n,1);
V_Minn = zeros(k*n,1);
sig2 = zeros(n,1);    
tmpY = [Y0(end-4+1:end,:); Yt];
U_hat = zeros(Tt,n);
for i=1:n
    Z = [ones(Tt,1) tmpY(4:end-1,i) tmpY(3:end-2,i) tmpY(2:end-3,i)...
        tmpY(1:end-4,i)];
    tmpb = (Z'*Z)\(Z'*tmpY(5:end,i));
    U_hat(:,i) = tmpY(5:end,i)-Z*tmpb;
    sig2(i) = mean(U_hat(:,i).^2);
end
Sig_hat = sig2;
count = 1;
for i=1:n
    for ii=0:k-1
        j = mod(ii,n); % variable index
        if j==0
            j = n;
        end
        l = ceil(ii/n); % lag length        
        if ii==0 % intercept
            V_Minn(count) = c3;
        elseif i==j % own lag
            V_Minn(count) = c1/l^2;
        elseif i~=j % lag of another variable        
            V_Minn(count) = c2*sig2(i)/(l^2*sig2(j));
        end
        count = count + 1;
    end
end

end