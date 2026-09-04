% bvar.priors.minn - Minnesota prior constructor for a VAR(p) with intercept.
% Extracted 2026-09-01 (step 4, SV/prior core). Canonicalizes:
%   chan2020_springer_largebvar/legacy/prior_Minn.m  -> call with n0pre = p
%   chan2023_joe_mlvarsv/legacy/utility/prior_Minn.m -> call with n0pre = 4
% Body is verbatim from the ml_varsv copy (the superset: it also returns the
% AR(4) residuals U_hat). The ONLY parameterization is n0pre, the number of
% presample rows of Y0 prepended before the univariate AR(4) fits:
%   tmpY = [Y0(end-n0pre+1:end,:); Yt]
% The two legacy copies differ only in that index (p vs hard-coded 4) and in
% the extra U_hat output; requesting three outputs reproduces the large_BVAR
% copy exactly (sig2 via stored residuals is bit-identical to the inline
% mean((y-Z*b).^2)). NOTE: the AR(4) design matrix is conformable only when
% the prepended block has exactly 4 rows, so n0pre = p runs only for p = 4
% (as in every large_BVAR caller); at p = 4 the two settings coincide.
%
% Prior: alpha ~ N(beta_Minn, diag(V_Minn)), equation by equation, with
% intercept variance c3, own-lag variance c1/l^2, cross-lag variance
% c2*sig2_i/(l^2*sig2_j); sig2 are univariate AR(4) residual variances.
%
% See:
% Chan, J.C.C. (2020). Large Bayesian Vector Autoregressions. In: P. Fuleky (Eds),
% Macroeconomic Forecasting in the Era of Big Data, 95-125, Springer, Cham
% Chan, J.C.C. (2023). Comparing Stochastic Volatility Specifications for
% Large Bayesian VARs, Journal of Econometrics, 235(2): 1419-1446.

function [beta_Minn,V_Minn,Sig_hat,U_hat] = minn(p,c1,c2,c3,Y0,Yt,n0pre)
[Tt,n] = size(Yt);
k = 1 + n*p;
beta_Minn = zeros(k*n,1);
V_Minn = zeros(k*n,1);
sig2 = zeros(n,1);
tmpY = [Y0(end-n0pre+1:end,:); Yt];
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
