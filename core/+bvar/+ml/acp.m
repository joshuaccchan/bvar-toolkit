% bvar.ml.acp - log marginal likelihood of a VAR under the asymmetric conjugate
% prior, in CLOSED FORM.
%
%   lml = bvar.ml.acp(p, Y, Z, prior)
%   lml = bvar.ml.acp(p, Y, Z, prior, 'ridge', 1e-6)
%
%   Z     : T x (n p + 1) lag matrix, intercept first (bvar.util.build_lags)
%   prior : the struct from bvar.priors.acp_redu or acp_stru
%   ridge : added to the diagonal of the posterior precision; default 0
%
% The ridge moves the value, so it is a modelling choice. By default the
% posterior precision is formed as iVi + Xi'*Xi; a positive ridge adds
% ridge*speye(ki) to it, a jitter that keeps the Cholesky alive at n = 35. On
% the ACP package's 15-variable dataset at its own kappa = (.04, .0016, 1, 100),
% ridge = 0 and ridge = 1e-6 give log marginal likelihoods 1.99 apart, which is
% large enough to affect a model comparison. Hold the setting fixed across the
% models being compared, and report which one was used.
%
% This is the property that motivates the prior. Conjugacy makes each equation's
% marginal likelihood a ratio of normal-inverse-gamma normalizing constants, so
% the whole quantity is a sum over equations of terms in log|Vi|, the Cholesky
% factor of the posterior precision, and gammaln - no simulation, no importance
% sampling, no second pass over stored draws. Selecting the lag length or the
% shrinkage hyperparameters by maximizing it is then a small optimization
% problem; bvar.priors.acp_opt_kappa does exactly that.
%
% Contrast the other entries in this namespace, which exist because their models
% have no such expression: the kron_bvar family needs Chib's method and the
% mlvarsv family needs adaptive importance sampling, both requiring the chain.
%
% Provenance and the legacy copies this stands in for: tests/variant_map.md.
% Equivalence: tests/unit/test_acp_equivalence.m.
%
% See:
% Chan, J.C.C. (2022). Asymmetric Conjugate Priors for Large Bayesian VARs,
% Quantitative Economics, 13(3): 1145-1169.

function lml = acp(p,Y,Z,prior,varargin)
ridge = 0;
for iv = 1:2:numel(varargin)
    switch lower(varargin{iv})
        case 'ridge', ridge = varargin{iv+1};
        otherwise, error('bvar:ml:acp:badOption', ...
                'unknown option ''%s''', varargin{iv});
    end
end
[T,n] = size(Y);
lml = -n*T/2*log(2*pi);
count_alp = 0;
for ii = 1:n
    yi = Y(:,ii);
    ki = n*p+ii;
    mi = [prior.beta0(:,ii);prior.alp0(count_alp+1:count_alp+ii-1)];
    Vi = sparse(1:ki,1:ki,[prior.Vbeta(:,ii);prior.Valp(count_alp+1:count_alp+ii-1)]);
    nui = prior.nu(ii);
    Si = prior.S(ii);
    Xi = [Z -Y(:,1:ii-1)];

    iVi = Vi\speye(ki);
    Kthetai = iVi + Xi'*Xi;
    if ridge ~= 0                       % branch, not + 0*speye: the default
        Kthetai = Kthetai + ridge*speye(ki);   % path's arithmetic is untouched
    end
    CKthetai = chol(Kthetai,'lower');
    thetai_hat = CKthetai'\(CKthetai\(iVi*mi + Xi'*yi));
    Si_hat = Si + (yi'*yi + mi'*iVi*mi - thetai_hat'*Kthetai*thetai_hat)/2;

    lml = lml -1/2*(sum(log(diag(Vi))) + 2*sum(log(diag(CKthetai)))) ...
        + nui*log(Si) - (nui+T/2)*log(Si_hat) + gammaln(nui+T/2) - gammaln(nui);

    count_alp = count_alp + ii-1;
end
end
