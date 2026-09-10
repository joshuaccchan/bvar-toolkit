% bvar.priors.acp_stru - asymmetric conjugate prior elicited directly on the
% structural parameterization of the VAR.
%
%   prior = bvar.priors.acp_stru(n, p, kappa, sig2)
%   prior = bvar.priors.acp_stru(n, p, kappa, sig2, idx_ns)
%
%   n, p   : number of variables and lag length
%   kappa  : [kappa1 kappa2 kappa3 kappa4] - own lags, other lags, impact
%            matrix, intercepts
%   sig2   : n-vector of AR(4) residual variances (bvar.priors.resid_var_ar4)
%   idx_ns : indices of nonstationary variables, whose first own lag gets prior
%            mean one (default none)
%   prior  : struct with fields beta0, Vbeta (VAR coefficients, (n p + 1) x n,
%            equation by equation), alp0, Valp (free elements of the impact
%            matrix, stacked equation by equation), nu, S (the n inverse-gamma
%            error-variance parameters)
%
% Do not merge the local helper prior_ACPi with the 2019 working paper copy of
% the same name: that one has no nonstationary unit-mean option on the first
% own lag, and returns Vi as a sparse diagonal matrix where this one returns a
% vector.
%
% Provenance and the legacy copies this stands in for: tests/variant_map.md.
% Equivalence: tests/unit/test_acp_stru.m, test_acp_2019wp_relationship.m.
%
% See:
% Chan, J.C.C. (2022). Asymmetric Conjugate Priors for Large Bayesian VARs,
% Quantitative Economics, 13(3): 1145-1169

function prior = acp_stru(n,p,kappa,sig2,idx_ns)
if nargin == 4
    idx_ns = [];
end
k_beta = n*(n*p+1);
k_alp = n*(n-1)/2;
prior.beta0 = zeros(k_beta/n,n);
prior.alp0 = zeros(k_alp,1);
prior.Vbeta = zeros(k_beta/n,n);
prior.Valp = zeros(k_alp,1);
prior.nu = zeros(n,1);
prior.S = zeros(n,1);
count_alp = 0;
for ii = 1:n
    is_ns = any(idx_ns == ii);
    [mi,Vi,nui,Si] = prior_ACPi(n,p,ii,kappa,sig2,is_ns);
    prior.beta0(:,ii) = mi(1:k_beta/n);
    prior.alp0(count_alp+1:count_alp+ii-1) = mi(k_beta/n+1:end);
    prior.Vbeta(:,ii) = Vi(1:k_beta/n);
    prior.Valp(count_alp+1:count_alp+ii-1) = Vi(k_beta/n+1:end);
    prior.nu(ii) = nui;
    prior.S(ii) = Si;
    count_alp = count_alp + ii - 1;
end

end

function [mi,Vi,nui,Si] = prior_ACPi(n,p,var_i,kappa,sig2,is_ns)
ki = var_i + n*p;
mi = zeros(ki,1);
Vi = zeros(ki,1);
    % construct Vi
for j=1:ki
    if j <= n*p+1
        l = ceil((j-1)/n); % lag length
        idx = mod(j-1,n);  % variable index
        if idx==0
            idx = n;
        end
    else
        idx = j - (n*p+1);
    end

    if j==1 % intercept
        Vi(j) = kappa(4);
    elseif j > n*p+1    % alpha_i
        Vi(j) = kappa(3)/sig2(idx);
    elseif idx == var_i % own lag
        Vi(j) = kappa(1)/(l^2*sig2(idx));
        if l == 1 && is_ns % if first own lag & variable is nonstationary
            mi(j) = 1;
        end
    else % lag of other variables
        Vi(j) = kappa(2)/(l^2*sig2(idx));
    end
end
Si = sig2(var_i)/2;
nui = 1 + var_i/2;
end
