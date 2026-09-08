% bvar.priors.minnesota_C - Minnesota-prior second-moment building blocks C and
% the index sets of the own-lag (kappa1) and other-lag (kappa2) coefficients,
% equation by equation (intercept first, then lag-1 block, lag-2 block, ...).
%
% Body from chan2021_ijf_mahp/legacy/get_C.m, renamed; the
% hybtvp, mlvarsv and oisv packages carry copies identical modulo comments.
% Equivalence: tests/unit/test_minnesota_C.m. Record: tests/variant_map.md.
%
% This function constructs the second moments of the Minnesota prior
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
