% bvar.samplers.gig_shrinkage - the kappa/psi generalized-inverse-Gaussian
% hierarchical shrinkage block of Chan (2021, IJF): draws the global shrinkage
% hyperparameters via gigrnd and, for the normal-gamma variants, the
% coefficient-specific local scales psi_kappa1/psi_kappa2 (one gigrnd call
% each, floored at psi_floor).
%
% The function has three explicitly NAMED variants; the blocks are numerically
% DIFFERENT across models - never unify them:
%   'mng'  -> kappa(1:2) with the Minnesota C, then the psi block;
%   'ng'   -> a single kappa, no Minnesota C and no factor 2 - the NG prior
%             variance is kappa*psi;
%   'minn' -> kappa draws only, no psi block: psi_kappa1/psi_kappa2 pass
%             through untouched, callers may pass [], and nu_psi/psi_floor are
%             not referenced.
% Caller contract: the Psi reassembly, Psi(idx_kappa1) = psi_kappa1 and
% Psi(idx_kappa2) = psi_kappa2, stays with the caller.
% NEVER-MERGE: the NG forecasting sampler forecast_BVAR_NG.m is NOT reproduced
% by 'ng' at any psi_floor - its conditionals carry an extra factor 2, pairing
% with its doubled Valp/Vbeta; functionize it separately if the forecast
% pipeline is ever consolidated.
%
% rng consumption (all draws through gigrnd, resolved from third_party/):
%   'mng' : 2 + n*p + (n-1)*n*p gigrnd calls, in that order;
%   'ng'  : 1 + n*p + (n-1)*n*p gigrnd calls;
%   'minn': 2 gigrnd calls.
%
% Inputs:  variant     - 'mng' | 'ng' | 'minn'
%          beta        - n*(n*p+1) x 1 current coefficient draw
%          idx_kappa1  - n*p x 1 own-lag positions (from bvar.priors.minnesota_C)
%          idx_kappa2  - (n-1)*n*p x 1 other-lag positions
%          C           - Minnesota second-moment vector ('mng'/'minn'; unused by 'ng')
%          kappa       - current kappa: 4-vector, elements 1:2 updated ('mng'/'minn');
%                        scalar ('ng')
%          psi_kappa1  - n*p x 1 local scales ('mng'/'ng'; pass [] for 'minn')
%          psi_kappa2  - (n-1)*n*p x 1 local scales (idem)
%          nu_psi      - normal-gamma shape ('mng'/'ng'; unused by 'minn')
%          c01, c02    - gamma prior [shape, rate] pairs (c02 unused by 'ng')
%          n, p        - VAR dimensions
%          psi_floor   - psi lower bound against arithmetic underflow
%                        ('mng'/'ng'; unused by 'minn'). The estimation
%                        samplers use 1e-10, the MNG forecasting sampler 1e-16
% Outputs: kappa, psi_kappa1, psi_kappa2 - updated state (psi pass through 'minn')
%
% See:
% Chan, J.C.C. (2021). Minnesota-Type Adaptive Hierarchical Priors for
% Large Bayesian VARs, International Journal of Forecasting, 37(3), 1212-1226.

function [kappa, psi_kappa1, psi_kappa2] = gig_shrinkage(variant, beta, ...
    idx_kappa1, idx_kappa2, C, kappa, psi_kappa1, psi_kappa2, nu_psi, ...
    c01, c02, n, p, psi_floor)
switch variant
    case 'mng'
            % sample kappa1 and kappa2
        tmpc1 = sum(beta(idx_kappa1).^2./(2*psi_kappa1.*C(idx_kappa1)));
        tmpc2 = sum(beta(idx_kappa2).^2./(2*psi_kappa2.*C(idx_kappa2)));
        kappa(1) = gigrnd(c01(1)-n*p/2,2*c01(2),tmpc1,1);
        kappa(2) = gigrnd(c02(1)-(n-1)*n*p/2,2*c02(2),tmpc2,1);

            % sample psi
        tmpv1 = beta(idx_kappa1).^2./(2*C(idx_kappa1)*kappa(1));
        tmpv2 = beta(idx_kappa2).^2./(2*C(idx_kappa2)*kappa(2));
        for ik=1:1:n*p  % lower bound psi_floor to avoid arithmetic underflow
            psi_kappa1(ik) = max(gigrnd(nu_psi-1/2,nu_psi,tmpv1(ik),1),psi_floor);
        end
        for il=1:1:(n-1)*n*p
            psi_kappa2(il) = max(gigrnd(nu_psi-1/2,nu_psi,tmpv2(il),1),psi_floor);
        end
    case 'ng'
            % sample kappa
        tmpc1 = sum(beta(idx_kappa1).^2./psi_kappa1);
        tmpc2 = sum(beta(idx_kappa2).^2./psi_kappa2);
        kappa = gigrnd(c01(1)-n^2*p/2,2*c01(2),tmpc1+tmpc2,1);

            % sample psi
        tmpv1 = beta(idx_kappa1).^2/kappa;
        tmpv2 = beta(idx_kappa2).^2/kappa;
        for ik=1:1:n*p  % lower bound psi_floor to avoid arithmetic underflow
            psi_kappa1(ik) = max(gigrnd(nu_psi-1/2,nu_psi,tmpv1(ik),1),psi_floor);
        end
        for il=1:1:(n-1)*n*p
            psi_kappa2(il) = max(gigrnd(nu_psi-1/2,nu_psi,tmpv2(il),1),psi_floor);
        end
    case 'minn'
            % sample kappa1 and kappa2
        tmpc1 = sum(beta(idx_kappa1).^2./C(idx_kappa1));
        tmpc2 = sum(beta(idx_kappa2).^2./C(idx_kappa2));
        kappa(1) = gigrnd(c01(1)-n*p/2,2*c01(2),tmpc1,1);
        kappa(2) = gigrnd(c02(1)-(n-1)*n*p/2,2*c02(2),tmpc2,1);
    otherwise
        error('bvar:samplers:gig_shrinkage', ...
            'unknown variant ''%s''; use mng, ng or minn', variant);
end
end
