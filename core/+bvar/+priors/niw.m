% bvar.priors.niw - natural-conjugate (normal-inverse-Wishart) prior constructor
% for a VAR(p) with intercept.
%
% Body from chan2023_joe_mlvarsv/legacy/utility/prior_NCP.m
% (the superset with U_hat), renamed, and parameterized by `variant` into four
% settings, each reproducing one legacy copy exactly:
%
%   'largebvar_nc'  -> chan2020_springer_largebvar/legacy/prior_NC.m
%       kappa = [c1 c2]; presample rows Y0(end-p+1:end,:); VA0 dense k x 1
%       vector; S0 = diag(sig2); nu0 = n+3.
%   'mlvarsv_ncp'   -> chan2023_joe_mlvarsv/legacy/utility/prior_NCP.m
%       kappa = [c1 c2]; presample rows Y0(end-4+1:end,:) (hard-coded 4);
%       VA0 dense vector; S0 = diag(sig2); nu0 = n+3; extra output U_hat.
%   'opthyper_ncp'  -> cjz2019_ad_opthyper/legacy/prior_NCP.m
%       kappa = [kappa1..kappa5]; presample rows Y0(end-p+1:end,:);
%       VA0(1) = kappa(3), VA0(i) = kappa(1)/(l^kappa(2)*sig2(idx));
%       S0 = kappa(5)*diag(sig2); nu0 = kappa(4)+n+1;
%       VA0 returned as SPARSE diagonal k x k matrix (as in the legacy copy).
%   'kron_script'   -> chan2020_jbes_kronecker/legacy/construct_prior_A.m
%       (a workspace SCRIPT; its callers all set S0 = eye(n); nu0 = n+3
%       immediately before running it, and the script hard-codes
%       c1 = .2^2, c2 = 100). Reproduced by kappa = [.2^2 100]:
%       A0/VA0/sig2 as the script computes them, plus S0 = eye(n), nu0 = n+3.
%
% NOTE: the AR(4) design matrix is conformable only when the prepended presample
% block has exactly 4 rows, so the Y0(end-p+1:end,:) variants run only for p = 4
% (as in all their legacy callers).
% Equivalence: tests/unit/test_prior_niw_<variant>.m, one per variant.
% Record: tests/variant_map.md.
%
% See:
% Chan, J.C.C. (2020). Large Bayesian Vector Autoregressions. In: P. Fuleky (Eds),
% Macroeconomic Forecasting in the Era of Big Data, 95-125, Springer, Cham
% Chan, J.C.C. (2023). Comparing Stochastic Volatility Specifications for
% Large Bayesian VARs, Journal of Econometrics, 235(2): 1419-1446.
% Chan, J. C. C., L. Jacobi, and D. Zhu (2020). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation,
% Journal of Forecasting, 39(6): 934-943.
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics,
% 38(1), 68-79.

function [A0,VA0,nu0,S0,U_hat] = niw(p,kappa,Y0,Y,variant)
switch variant
    case {'largebvar_nc','mlvarsv_ncp','opthyper_ncp','kron_script'}
    otherwise
        error('bvar:priors:niw:unknownVariant', ...
            'unknown variant ''%s''; use largebvar_nc, mlvarsv_ncp, opthyper_ncp or kron_script', variant);
end
[Tt,n] = size(Y);
k = 1+n*p;
A0 = zeros(k,n);
VA0 = zeros(k,1);
sig2 = zeros(n,1);
U_hat = zeros(Tt,n);
    % construct VA0
if strcmp(variant,'mlvarsv_ncp')
    tmpY = [Y0(end-4+1:end,:); Y];
else
    tmpY = [Y0(end-p+1:end,:); Y];
end
for i=1:n
    Z = [ones(Tt,1) tmpY(4:end-1,i) tmpY(3:end-2,i) tmpY(2:end-3,i)...
        tmpY(1:end-4,i)];
    tmpb = (Z'*Z)\(Z'*tmpY(5:end,i));
    U_hat(:,i) = tmpY(5:end,i)-Z*tmpb;
    sig2(i) = mean(U_hat(:,i).^2);
end
if strcmp(variant,'opthyper_ncp')
    for i=1:k
        l = ceil((i-1)/n);
        idx = mod(i-1,n); % variable index
        if idx==0
            idx = n;
        end
        if i==1 % intercept
            VA0(1) = kappa(3);
        else
            VA0(i) = kappa(1)/(l^kappa(2)*sig2(idx));
        end
    end
    S0 = kappa(5)*diag(sig2); nu0 = kappa(4)+n+1;
    VA0 = sparse(1:k,1:k,VA0);
else
    c1 = kappa(1); c2 = kappa(2);
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
    nu0 = n+3;
    if strcmp(variant,'kron_script')
        S0 = eye(n);
    else
        S0 = diag(sig2);
    end
end
end
