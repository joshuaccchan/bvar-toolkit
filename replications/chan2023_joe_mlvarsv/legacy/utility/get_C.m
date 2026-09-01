% This function constructs the second moments of the Minnesota prior
%
% See:
% Chan, J.C.C. (2022). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, forthcoming.

function [C,idx_kappa1,idx_kappa2] = get_C(n,p,sig2)
k_beta = n^2*p+n;
C = zeros(k_beta,1);
idx_kappa1 = [];
idx_kappa2 = [];
count = 1;

for ii = 1:n
    Ci = zeros(n*p+1,1);
        % construct Ci
    for j=1:n*p+1    
        l = ceil((j-1)/n); % lag length
        idx = mod(j-1,n);  % variable index
        if idx==0
            idx = n;
        end         
        if j==1 % intercept            
            Ci(j) = sig2(ii);
        elseif idx == ii % own lag
            Ci(j) = 1/l^2;
            idx_kappa1 = [idx_kappa1; count];
        else % lag of other variables            
            Ci(j) = sig2(ii)/(l^2*sig2(idx));
            idx_kappa2 = [idx_kappa2; count];
        end 
            count = count + 1;
    end
    C((ii-1)*(n*p+1)+1:ii*(n*p+1)) = Ci;
end
end