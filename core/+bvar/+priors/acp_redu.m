% bvar.priors.acp_redu - asymmetric conjugate prior elicited on the
% reduced-form parameterization, mapped to the implied prior on the
% structural parameterization.
% Extracted 2026-09-01 (step 4, SV/prior core). Canonical source:
% replications/chan2022_qe_acp/legacy/utility/prior_ACP_redu.m (R1, the
% published QE version). Body verbatim; only two edits: the function was
% renamed prior_ACP_redu -> acp_redu, and the internal call
% prior_ACP_stru(...) now targets the canonicalized core copy
% bvar.priors.acp_stru(...) (itself verbatim from the same R1 package).
% No parameterization was added.
%
% This function first elicits the asymmetric conjugate prior on the
% reduced-form parameterization and then constructs the implied prior on
% the structural parameterization
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
