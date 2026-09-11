% bvar.priors.minnesota_C - Minnesota-prior second-moment building blocks C and
% the index sets of the own-lag (kappa1) and other-lag (kappa2) coefficients,
% equation by equation (intercept first, then lag-1 block, lag-2 block, ...).
%
%   [C,idx_kappa1,idx_kappa2] = bvar.priors.minnesota_C(n, p, sig2)
%
%   n, p       : number of variables and lag length
%   sig2       : n-vector of AR(4) residual variances
%   C          : (n^2 p + n) x 1 vector of second moments - sig2_i on the
%                intercepts, 1/l^2 on own lags, sig2_i/(l^2*sig2_j) on other
%                lags. The shrinkage hyperparameters scale it in
%                bvar.priors.vtheta.
%   idx_kappa1 : indices into C of the own-lag coefficients
%   idx_kappa2 : indices into C of the other-lag coefficients
%
% See:
% Chan, J.C.C. (2021). Minnesota-Type Adaptive Hierarchical Priors for
% Large Bayesian VARs, International Journal of Forecasting, 37(3): 1212-1226

function [C,idx_kappa1,idx_kappa2] = minnesota_C(n,p,sig2)
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
