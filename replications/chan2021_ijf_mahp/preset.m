% chan2021_ijf_mahp/preset - every constant the legacy MAHP ESTIMATION pipeline
% hard-codes (main_BVAR.m dispatching to the workspace scripts BVAR_MNG.m /
% BVAR_NG.m / BVAR_Minn.m), one field per constant, each cited to its legacy
% source line. Consumed by run_all.m in this folder. Compiled 2026-09-01
% (step 5, MAHP flagship functionization); the legacy folder is never modified -
% this file only transcribes it.
%
% The FORECAST pipeline (main_forecasting.m + forecast_BVAR_*.m) is NOT
% functionized in this step; its known numerical divergences from the
% estimation scripts (2026-09-01 audit) are recorded under pr.forecast as
% documentation-only fields, so a future forecast functionization does not
% silently inherit the estimation settings (or vice versa).
%
% See:
% Chan, J.C.C. (2021). Minnesota-Type Adaptive Hierarchical Priors for
% Large Bayesian VARs, International Journal of Forecasting, 37(3), 1212-1226.

function pr = preset()

    % ---- data & model dimensions (main_BVAR.m) ----
pr.data_file  = 'macrodata_Q_2018Q4.csv';   % main_BVAR.m line 21; lives in legacy/, loaded read-only by run_all
pr.data_rows  = 238;                        % main_BVAR.m line 22: data = macrodata_Q_2018Q4(1:238,:)
pr.var_id     = [1,2,22,23,35,37,57,58,59,76,81,83,95,120,138,144,145,147,148,152,160,161,245]; % main_BVAR.m line 23 (n = 23)
pr.var_id_n30 = [1,2,22,23,35,37,57,58,59,76,81,83,95,120,138,144,145,147,148,152,...
    160,161,245,3,18,97,123,133,156,199];   % main_BVAR.m lines 24-25: n = 30 alternative, COMMENTED OUT in legacy (documented only; run_all uses pr.var_id)
pr.n0         = 8;                          % main_BVAR.m lines 26-27: first 8 obs are the initial conditions Y0
pr.p          = 4;                          % main_BVAR.m line 19 (legacy comment: if p > 8, Y0/Y must change)
pr.nsim_default   = 10000;                  % main_BVAR.m line 18
pr.burnin_default = 1000;                   % main_BVAR.m line 18

n = numel(pr.var_id);                       % = 23; legacy gets it as [T,n] = size(Y), main_BVAR.m line 29

    % ---- prior hyperparameters (main_BVAR.m lines 42-46) ----
pr.ah   = zeros(n,1);                       % line 42: prior mean of the initial log-volatilities h0
pr.Vh   = 10*ones(n,1);                     % line 42: prior variance of h0
pr.nuh0 = 5*ones(n,1);                      % line 43: IG shape for the SV state variances Sigh
pr.Sh0  = .01*ones(n,1).*(pr.nuh0-1);       % line 43: IG scale for Sigh (= .04); also the Sigh chain INIT (BVAR_MNG.m line 21, BVAR_NG.m line 20, BVAR_Minn.m line 12)
pr.c01  = [1, 1/.04];                       % line 44: gamma prior [shape, rate] for kappa1 (MNG/Minn) / the single kappa (NG)
pr.c02  = [1, 1/.04^2];                     % line 45: gamma prior [shape, rate] for kappa2 (1/.04^2 = 625)
pr.lam0_nu_psi = 1;                         % line 46: exponential prior rate for nu_psi (MNG/NG)

    % ---- constants shared by the three estimation samplers ----
pr.sv_offset = .0001;                       % Ystar = log(U.^2 + .0001): BVAR_MNG.m line 63, BVAR_NG.m line 61, BVAR_Minn.m line 54
pr.psi_floor = 1e-10;                       % psi lower bound (underflow guard): BVAR_MNG.m lines 77/80, BVAR_NG.m lines 74/77

    % ---- BVAR_MNG.m (model 1: Minnesota-type normal-gamma) ----
pr.mng.nu_psi_init  = .5;                   % BVAR_MNG.m line 17
pr.mng.kappa_init   = [.4,.001,1,100];      % BVAR_MNG.m line 18: [kappa1 own-lag, kappa2 other-lag, kappa3 impact, kappa4 intercept]
pr.mng.vtheta_scale = 2;                    % BVAR_MNG.m line 39: Valp = 2*Valp; Vbeta = 2*Vbeta (prior variance is 2*kappa*psi*C)
    % psi chain init (BVAR_MNG.m lines 24-25):
    %   psi_kappa1 ~ gamrnd(nu_psi, 2/nu_psi, n*p, 1)
    %   psi_kappa2 ~ gamrnd(nu_psi, 2/nu_psi, (n-1)*n*p, 1)
    % run_all keeps those expressions verbatim; the only free constant is
    % pr.mng.nu_psi_init above (scale 2/nu_psi = 4 at the init value).

    % ---- BVAR_NG.m (model 2: plain normal-gamma) ----
pr.ng.nu_psi_init = .5;                     % BVAR_NG.m line 16
pr.ng.kappa_init  = .04;                    % BVAR_NG.m line 17: single global kappa (scalar)
pr.ng.kappa3      = 1;                      % BVAR_NG.m line 37: getVtheta(...,[kappa,kappa,1,100],Psi,sig2) - impact-matrix scale
pr.ng.kappa4      = 100;                    % BVAR_NG.m line 37: intercept scale
    % NB (BVAR_NG.m line 37): estimation BVAR_NG does NOT rescale Valp/Vbeta
    % (contrast pr.mng.vtheta_scale, and pr.forecast.ng_doubles_vtheta below),
    % and it passes Psi - not C.*Psi - to getVtheta: the NG prior variance is
    % kappa*psi with no Minnesota C scaling.
    % psi chain init (BVAR_NG.m lines 23-24): identical expressions to BVAR_MNG.

    % ---- BVAR_Minn.m (model 3: data-based Minnesota) ----
pr.minn.kappa_init = [.04,.04,1,100];       % BVAR_Minn.m line 9 - NOT the MNG init [.4,.001,1,100]
    % No psi block and no nu_psi step in this model.

    % ---- estimation-vs-forecast divergences (2026-09-01 audit) ----
    % DOCUMENTATION ONLY: the forecast pipeline is not functionized in this
    % step; run_all reads nothing below this line.
pr.forecast.nsim_default   = 20000;         % main_forecasting.m line 19
pr.forecast.burnin_default = 100;           % main_forecasting.m line 20 (vs estimation 1000)
pr.forecast.T0 = 91;                        % main_forecasting.m line 34: first forecast origin (legacy comment: 1974Q1)
pr.forecast.c03 = [1, 1/1];                 % main_forecasting.m line 41 - defined but UNUSED by all three forecast samplers
pr.forecast.c04 = [1, 1/100];               % main_forecasting.m line 42 - defined but UNUSED
pr.forecast.psi_floor = 1e-16;              % forecast_BVAR_MNG.m lines 79/82, forecast_BVAR_NG.m lines 80/83 - NOT the estimation 1e-10
pr.forecast.psi1_init_scale_halved = true;  % forecast_BVAR_MNG.m line 31, forecast_BVAR_NG.m line 32:
                                            %   psi_kappa1 ~ gamrnd(nu_psi, 2/nu_psi/2, n*p, 1) - scale HALVED vs
                                            %   estimation; psi_kappa2 keeps scale 2/nu_psi in both pipelines
pr.forecast.minn_kappa_init = [.4,.001,1,100]; % forecast_BVAR_Minn.m line 27 - the MNG init, NOT the estimation
                                            %   Minn init [.04,.04,1,100] (MNG/NG kappa inits agree across pipelines)
pr.forecast.ng_doubles_vtheta = true;       % forecast_BVAR_NG.m line 43: Valp = 2*Valp; Vbeta = 2*Vbeta - ABSENT from
                                            %   estimation BVAR_NG (forecast_BVAR_MNG line 42 doubles like estimation MNG line 39)
pr.forecast.ng_conditionals_carry_half = true; % forecast_BVAR_NG.m lines 72-78: tmpc_j = sum(beta_j.^2./(2*psi_j)),
                                            %   tmpv_j = beta_j.^2/(2*kappa) - the factor 2 pairs with the doubled
                                            %   Vtheta; estimation BVAR_NG has no factor 2. NEVER-MERGE:
                                            %   bvar.samplers.gig_shrinkage('ng',...) reproduces the ESTIMATION block only.
    % Further forecast-only mechanics with no draw-sequence divergence beyond
    % the above: sig2/C are recomputed from each expanding vintage (Y0,Yt);
    % only h(end,:) is stored; forecast_BVAR_MNG calls sample_nu_psi with one
    % output while forecast_BVAR_NG captures [nu_psi, flag] but never
    % accumulates the flag.
end
