% bvt.samplers.gig_shrinkage - the kappa/psi generalized-inverse-Gaussian
% shrinkage-update ladder of Chan (2021, IJF): draws the global shrinkage
% hyperparameters via gigrnd and, for the normal-gamma variants, the
% coefficient-specific local scales psi_kappa1/psi_kappa2 (one gigrnd call
% each, floored at psi_floor).
%
% Extracted 2026-09-01 (step 5, MAHP flagship functionization). One function,
% three explicitly NAMED variants; each branch is verbatim from its legacy
% source (the blocks are numerically DIFFERENT across models - never unify):
%   variant 'mng'  -> chan2021_ijf_mahp/legacy/BVAR_MNG.m  lines 68-81:
%                     kappa(1) ~ GIG(c01(1)-n*p/2,       2*c01(2), sum(beta1.^2./(2*psi1.*C1)))
%                     kappa(2) ~ GIG(c02(1)-(n-1)*n*p/2, 2*c02(2), sum(beta2.^2./(2*psi2.*C2)))
%                     psi_j(i) = max(GIG(nu_psi-1/2, nu_psi, beta_j(i)^2/(2*C_j(i)*kappa_j)), psi_floor)
%   variant 'ng'   -> chan2021_ijf_mahp/legacy/BVAR_NG.m   lines 66-78:
%                     single kappa ~ GIG(c01(1)-n^2*p/2, 2*c01(2), sum(beta1.^2./psi1)+sum(beta2.^2./psi2))
%                     psi_j(i) = max(GIG(nu_psi-1/2, nu_psi, beta_j(i)^2/kappa), psi_floor)
%                     (no Minnesota C, no factor 2 - the NG prior variance is kappa*psi)
%   variant 'minn' -> chan2021_ijf_mahp/legacy/BVAR_Minn.m lines 59-62:
%                     kappa(1)/kappa(2) draws only, chi = sum(beta_j.^2./C_j);
%                     no psi ladder (psi_kappa1/psi_kappa2 pass through untouched;
%                     callers may pass [], and nu_psi/psi_floor are not referenced).
% Documented settings reproducing each legacy copy exactly:
%   estimation BVAR_MNG : 'mng' with psi_floor = 1e-10 (its lines 77/80);
%   estimation BVAR_NG  : 'ng'  with psi_floor = 1e-10 (its lines 74/77);
%   estimation BVAR_Minn: 'minn' (floor unused);
%   forecast_BVAR_MNG   : 'mng' with psi_floor = 1e-16 (its lines 79/82; the
%                         kappa/psi conditionals are otherwise identical);
%   forecast_BVAR_Minn  : 'minn' (identical block, its lines 65-68).
% NEVER-MERGE: forecast_BVAR_NG.m is NOT reproduced by 'ng' at any psi_floor -
% its conditionals carry an extra factor 2 (tmpc_j = sum(beta_j.^2./(2*psi_j)),
% tmpv_j = beta_j.^2/(2*kappa), lines 72-78), pairing with its doubled
% Valp/Vbeta (line 43); a numerically different parameterization, to be
% functionized separately if the forecast pipeline is ever consolidated.
% Edits made: wrapped as a function; kappa/psi state passed in and returned;
% the hard-coded 1e-10 floor promoted to the argument psi_floor (settings
% above); the Psi(idx) reassembly (BVAR_MNG lines 82-83) stays with the caller.
% Branch bodies otherwise verbatim, including the 1:1:n*p loop stride.
%
% rng consumption (all draws through gigrnd, resolved from third_party/):
%   'mng' : 2 + n*p + (n-1)*n*p gigrnd calls, in that order;
%   'ng'  : 1 + n*p + (n-1)*n*p gigrnd calls;
%   'minn': 2 gigrnd calls.
%
% Inputs:  variant     - 'mng' | 'ng' | 'minn'
%          beta        - n*(n*p+1) x 1 current coefficient draw
%          idx_kappa1  - n*p x 1 own-lag positions (from bvt.priors.minnesota_C)
%          idx_kappa2  - (n-1)*n*p x 1 other-lag positions
%          C           - Minnesota second-moment vector ('mng'/'minn'; unused by 'ng')
%          kappa       - current kappa: 4-vector, elements 1:2 updated ('mng'/'minn');
%                        scalar ('ng')
%          psi_kappa1  - n*p x 1 local scales ('mng'/'ng'; pass [] for 'minn')
%          psi_kappa2  - (n-1)*n*p x 1 local scales (idem)
%          nu_psi      - normal-gamma shape ('mng'/'ng'; unused by 'minn')
%          c01, c02    - gamma prior [shape, rate] pairs (c02 unused by 'ng')
%          n, p        - VAR dimensions
%          psi_floor   - psi lower bound ('mng'/'ng'; unused by 'minn')
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
            % sample kappa1 and kappa2    [BVAR_MNG.m lines 68-71]
        tmpc1 = sum(beta(idx_kappa1).^2./(2*psi_kappa1.*C(idx_kappa1)));
        tmpc2 = sum(beta(idx_kappa2).^2./(2*psi_kappa2.*C(idx_kappa2)));
        kappa(1) = gigrnd(c01(1)-n*p/2,2*c01(2),tmpc1,1);
        kappa(2) = gigrnd(c02(1)-(n-1)*n*p/2,2*c02(2),tmpc2,1);

            % sample psi    [BVAR_MNG.m lines 74-81]
        tmpv1 = beta(idx_kappa1).^2./(2*C(idx_kappa1)*kappa(1));
        tmpv2 = beta(idx_kappa2).^2./(2*C(idx_kappa2)*kappa(2));
        for ik=1:1:n*p  % lower bound psi_floor to avoid arithmetic underflow
            psi_kappa1(ik) = max(gigrnd(nu_psi-1/2,nu_psi,tmpv1(ik),1),psi_floor);
        end
        for il=1:1:(n-1)*n*p
            psi_kappa2(il) = max(gigrnd(nu_psi-1/2,nu_psi,tmpv2(il),1),psi_floor);
        end
    case 'ng'
            % sample kappa    [BVAR_NG.m lines 66-68]
        tmpc1 = sum(beta(idx_kappa1).^2./psi_kappa1);
        tmpc2 = sum(beta(idx_kappa2).^2./psi_kappa2);
        kappa = gigrnd(c01(1)-n^2*p/2,2*c01(2),tmpc1+tmpc2,1);

            % sample psi    [BVAR_NG.m lines 71-78]
        tmpv1 = beta(idx_kappa1).^2/kappa;
        tmpv2 = beta(idx_kappa2).^2/kappa;
        for ik=1:1:n*p  % lower bound psi_floor to avoid arithmetic underflow
            psi_kappa1(ik) = max(gigrnd(nu_psi-1/2,nu_psi,tmpv1(ik),1),psi_floor);
        end
        for il=1:1:(n-1)*n*p
            psi_kappa2(il) = max(gigrnd(nu_psi-1/2,nu_psi,tmpv2(il),1),psi_floor);
        end
    case 'minn'
            % sample kappa1 and kappa2    [BVAR_Minn.m lines 59-62]
        tmpc1 = sum(beta(idx_kappa1).^2./C(idx_kappa1));
        tmpc2 = sum(beta(idx_kappa2).^2./C(idx_kappa2));
        kappa(1) = gigrnd(c01(1)-n*p/2,2*c01(2),tmpc1,1);
        kappa(2) = gigrnd(c02(1)-(n-1)*n*p/2,2*c02(2),tmpc2,1);
    otherwise
        error('bvt:samplers:gig_shrinkage', ...
            'unknown variant ''%s''; use mng, ng or minn', variant);
end
end
