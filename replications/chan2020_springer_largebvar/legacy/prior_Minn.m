% This script constructs the Minnesota prior
%
% See:
% Chan, J.C.C. (2020). Large Bayesian Vector Autoregressions. In: P. Fuleky (Eds),
% Macroeconomic Forecasting in the Era of Big Data, 95-125, Springer, Cham

function [beta_Minn,V_Minn,Sig_hat] = prior_Minn(p,c1,c2,c3,Y0,shortYt)
[Tt,n] = size(shortYt);
k = 1+ n*p;
beta_Minn = zeros(k*n,1);
V_Minn = zeros(k*n,1);
sig2 = zeros(n,1);    
tmpY = [Y0(end-p+1:end,:); shortYt];
for i=1:n
    Z = [ones(Tt,1) tmpY(4:end-1,i) tmpY(3:end-2,i) tmpY(2:end-3,i)...
        tmpY(1:end-4,i)];
    tmpb = (Z'*Z)\(Z'*tmpY(5:end,i));
    sig2(i) = mean((tmpY(5:end,i)-Z*tmpb).^2);
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