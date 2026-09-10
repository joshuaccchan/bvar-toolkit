% bvar.priors.acp_redu - asymmetric conjugate prior elicited on the
% reduced-form parameterization, mapped to the implied prior on the
% structural parameterization.
%
%   prior = bvar.priors.acp_redu(n, p, kappa, sig2)
%   prior = bvar.priors.acp_redu(n, p, kappa, sig2, idx_ns)
%
%   n, p   : number of variables and lag length
%   kappa  : [kappa1 kappa2 kappa3 kappa4] - own lags, other lags, impact
%            matrix, intercepts
%   sig2   : n-vector of AR(4) residual variances (bvar.priors.resid_var_ar4)
%   idx_ns : indices of nonstationary variables, whose first own lag gets prior
%            mean one (default none)
%   prior  : struct with fields beta0, Vbeta (VAR coefficients, (n p + 1) x n,
%            equation by equation), alp0, Valp (free elements of the impact
%            matrix), nu, S (the n inverse-gamma error-variance parameters).
%            Only Vbeta differs from bvar.priors.acp_stru; the other fields are
%            passed through unchanged.
%
% Provenance and the legacy copies this stands in for: tests/variant_map.md.
% Equivalence: tests/unit/test_acp_redu.m.
%
% See:
% Chan, J.C.C. (2022). Asymmetric Conjugate Priors for Large Bayesian VARs,
% Quantitative Economics, 13(3): 1145-1169

function prior = acp_redu(n,p,kappa,sig2,idx_ns)
if nargin == 4
    idx_ns = [];
end
k_beta = n*(n*p+1);
prior_stru = bvar.priors.acp_stru(n,p,kappa,sig2,idx_ns);
prior.alp0 = prior_stru.alp0;
prior.beta0 = prior_stru.beta0;
prior.Valp = prior_stru.Valp;
prior.nu = prior_stru.nu;
prior.S = prior_stru.S;
prior.Vbeta = zeros(k_beta/n,n);
for ii=1:n
    for jj=1:n*p+1
        if ii == 1
            prior.Vbeta(jj,ii) = prior_stru.Vbeta(jj,ii);
        else
            prior.Vbeta(jj,ii) = prior_stru.Vbeta(jj,ii) ...
                + sum(prior_stru.Vbeta(jj,1:ii-1) ...
                + prior_stru.beta0(jj,1:ii-1).^2./sig2(1:ii-1)');
        end
    end
end

end
