% This function first elicits the asymmetric conjugate prior on the
% reduced-form parameterization and then constructs the implied prior on
% the structural parameterization
%
% See:
% Chan, J.C.C. (2021). Asymmetric conjugate priors for large Bayesian VARs,
% Quantitative Economics, forthcoming.

function prior = prior_ACP_redu(n,p,kappa,sig2,idx_ns)
if nargin == 4
    idx_ns = [];
end
k_beta = n*(n*p+1);
prior_stru = prior_ACP_stru(n,p,kappa,sig2,idx_ns);
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