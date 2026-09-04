% bvar.samplers.nu_psi_ng - independence-MH draw of the normal-gamma shape
% parameter nu_psi in the (Minnesota-type) normal-gamma priors of Chan (2021,
% IJF), given the local scales psi_kappa1/psi_kappa2 (psi ~ G(nu_psi, 2/nu_psi))
% under an Exp(lam0_nu_psi) prior: Newton-Raphson search for the conditional
% mode nu_psit (fminbnd fallback if an iterate goes negative), Gaussian
% N(nu_psit, -1/H) proposal, MH accept-reject.
% rng consumption: one randn always, then one rand IFF the candidate is
% positive (Newton/fminbnd are deterministic).
%
% Extracted 2026-09-01 (step 5, MAHP flagship functionization). Canonical
% source: chan2021_ijf_mahp/legacy/sample_nu_psi.m - the ONLY copy in the
% package; it already is a function (not an inline workspace block). Called by
% BVAR_MNG.m line 86, BVAR_NG.m line 83, and forecast_BVAR_NG.m line 89 with
% two outputs (the forecast call captures flag but never accumulates it), and
% by forecast_BVAR_MNG.m line 88 with one output - all four call sites are
% served by this single function unchanged.
% Edits made: provenance header prepended; function renamed
% sample_nu_psi -> nu_psi_ng. Body verbatim, including the inert `count`
% logic: count is initialized to 0 and never incremented, so the
% `count < 100` condition in the while loop and the `if count == 100` fminbnd
% fallback are dead code - kept verbatim rather than cleaned up.
%
% Inputs:  psi_kappa1, psi_kappa2 - local scale draws (n*p x 1, (n-1)*n*p x 1)
%          nu_psi                 - current draw
%          lam0_nu_psi            - exponential prior rate
% Outputs: nu_psi, flag (1 if accepted), f_nu_psi (handle: log target kernel)
%
% This script samples the parameter nu_psi in the Minnesota-type
% normal-gamma prior
%
% See:
% Chan, J.C.C. (2021). Minnesota-Type Adaptive Hierarchical Priors for
% Large Bayesian VARs, International Journal of Forecasting, 37(3), 1212-1226.

function [nu_psi, flag, f_nu_psi] = nu_psi_ng(psi_kappa1,psi_kappa2,...
    nu_psi,lam0_nu_psi)
flag = 0;
n2p = size(psi_kappa1,1) + size(psi_kappa2,1);
sum1 = sum(log(psi_kappa1)) + sum(log(psi_kappa2));
sum2 = sum(psi_kappa1) + sum(psi_kappa2);
f_nu_psi = @(x) n2p*(x*log(x/2) - gammaln(x)) + (x-1)*sum1 - x/2*sum2 - lam0_nu_psi*x;
df_nu_psi = @(x) n2p*(log(x/2) + 1 - psi(x)) + sum1 - .5*sum2 - lam0_nu_psi;
d2f_nu_psi = @(x) n2p/x - n2p*psi(1,x);
S_nu_psi = 1;
nu_psit = nu_psi;
count = 0;
while abs(S_nu_psi) > 1e-5 && count < 100  % stopping criteria
    S_nu_psi = df_nu_psi(nu_psit);
    H_nu_psi = d2f_nu_psi(nu_psit);
    nu_psit = nu_psit - H_nu_psi\S_nu_psi;
    if nu_psit < 0
        nu_psit = fminbnd(@(x)-f_nu_psi(x),.01,5);
        H_nu_psi = d2f_nu_psi(nu_psit);
        break;
    end
end
if count == 100
    nu_psit = fminbnd(@(x)-f_nu_psi(x),.01,5);
    H_nu_psi = d2f_nu_psi(nu_psit);
end

Dnu_psi = -1/H_nu_psi;
nu_psic = nu_psit + sqrt(Dnu_psi)*randn;
if nu_psic > 0
    alp_MH = f_nu_psi(nu_psic) - f_nu_psi(nu_psi) ...
        - .5*(nu_psi-nu_psit)^2/Dnu_psi + .5*(nu_psic-nu_psit)^2/Dnu_psi;
    if exp(alp_MH) > rand
        nu_psi = nu_psic;
        flag = 1;
    end
end
end
