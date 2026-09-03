% chan2020_jbes_kronecker/preset - every constant the legacy Kronecker
% FULL-SAMPLE ESTIMATION + MARGINAL-LIKELIHOOD pipeline hard-codes
% (main_BVAR.m dispatching to the workspace scripts BVAR*.m, each of which
% chains into its ml_BVAR_*.m when cp_ml = 1), one field per constant, each
% cited to its legacy source line. Consumed by run_all.m and run_ml.m in
% this folder. Compiled 2026-09-02 (step 8, Kronecker family pass, part 1);
% the legacy folder is never modified - this file only transcribes it. The
% REALTIME FORECASTING pipeline (main_forecasting.m + realtime_forecasts/)
% is part 2 and is NOT functionized in this step; only its path hazard is
% recorded below because it endangers the ESTIMATION/ML path too.
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics,
% 38(1), 68-79.

function pr = preset()

    % ---- run configuration (main_BVAR.m) ----
pr.model_default  = 8;                  % main_BVAR.m line 19 (BVAR-CSV-t-MA)
pr.cp_ml_default  = 1;                  % line 20: compute the marginal likelihood after estimation
pr.p              = 4;                  % line 21 (legacy comment: if p > 4, Y0 and shortY must change)
pr.nsims_default  = 30000;              % line 22
pr.burnin_default = 5000;               % line 23

    % ---- data (main_BVAR.m lines 26-31) ----
pr.data_file = 'data_Q.csv';            % line 26; 1959Q1-2013Q4, lives in legacy/, loaded read-only by run_all
pr.var_cols  = [1:3 6:15 17 19:24];     % line 27: the 20 variables of Appendix B out of the file's 24 columns
pr.n0        = 4;                       % line 28: first 4 obs are the initial conditions Y0

    % ---- VAR coefficient prior (construct_prior_A.m; S0/nu0 from the
    %      callers' prior blocks) ----
pr.minn_kappa = [.2^2 100];             % construct_prior_A.m line 12: c1 = .2^2 (slopes), c2 = 100 (intercepts)
                                        %   -> bvt.priors.niw(p, pr.minn_kappa, Y0, shortY, 'kron_script')
                                        %   S0 = eye(n), nu0 = n+3: set by every BVAR*.m immediately before
                                        %   running construct_prior_A (e.g. BVAR.m lines 9-10, BVAR_CSV.m line 9);
                                        %   the niw 'kron_script' variant returns exactly those.

    % ---- shared MCMC hyperparameters (per-script prior blocks) ----
pr.nuub = 100;                          % nu upper bound, t models: BVAR_t.m line 11, BVAR_t_CSV.m line 12,
                                        %   BVAR_t_MA.m line 12, BVAR_CSV_t_MA.m line 14
pr.nuh0 = 5;                            % IG shape for sigh2: BVAR_CSV.m line 11, BVAR_t_CSV.m line 10,
                                        %   BVAR_CSV_MA.m line 12, BVAR_CSV_t_MA.m line 12
pr.Sh0  = .01*(pr.nuh0-1);              % IG scale for sigh2 (= .04), same lines
pr.rho0 = .9;                           % rho prior mean, same CSV scripts (next line in each)
pr.Vrho = .2^2;                         % rho prior variance, same lines
pr.psi0 = 0;                            % psi prior mean, MA models: BVAR_MA.m line 9, BVAR_t_MA.m line 9,
                                        %   BVAR_CSV_MA.m line 9, BVAR_CSV_t_MA.m line 9
pr.Vpsi = 1;                            % psi prior variance, same lines; the lpri_psi handle adds a -1e10
                                        %   penalty outside (-.99, .99) (the following line in each script)
pr.progress_every = 5000;               % loop-counter disp cadence, every MCMC script

    % ---- chain initial values (per-script "initialize the chain" blocks) ----
    % h = zeros(T,1); rho = .8; sigh2 = .1; nu = 5; lam ~ 1./gamrnd(nu/2,2/nu,T,1);
    % psi = .1 (BVAR_MA.m line 31) but psi1 = -.1 (BVAR_t_MA.m line 40) and
    % psi = -.1 (BVAR_CSV_MA.m line 36, BVAR_CSV_t_MA.m line 41); psihat
    % initialized at the same value. run_all keeps these expressions verbatim
    % per model; they carry no free constants beyond the values above.

    % ---- rho MH truncation divergences (2026-09-02 audit; all VERIFIED
    %      against source, values kept verbatim in run_all / bvt.ml.*) ----
pr.rho_mh_bnd_est = [NaN NaN .9999 NaN .9999 NaN .9999 .99];
                                        % ESTIMATION scripts, by model number:
                                        %   .9999 in BVAR_CSV.m line 69, BVAR_t_CSV.m line 87,
                                        %   BVAR_CSV_MA.m line 83; but .99 in BVAR_CSV_t_MA.m line 92 -
                                        %   model 8's estimation bound is TIGHTER than its siblings'.
pr.rho_mh_bnd_ml  = [NaN NaN .9999 NaN .999 NaN .9999 .999];
                                        % ML reduced runs, by model number:
                                        %   ml_BVAR_CSV.m line 95 = .9999 (matches its estimation);
                                        %   ml_BVAR_t_CSV.m line 81 = .999 (estimation .9999);
                                        %   ml_BVAR_CSV_MA.m line 74 = .9999 (matches);
                                        %   ml_BVAR_CSV_t_MA.m line 98 = .999 (estimation .99).
                                        %   All four rho ordinate grids span linspace(-.999,.999,...)
                                        %   regardless of the MH bound in force.
pr.psi_mh_bnd = .99;                    % psi MH acceptance bound |psic| < .99 and fminbnd range
                                        %   (-.99,.99) in ALL MA scripts, estimation and ml alike;
                                        %   Kpsic fallback 1/.05^2 when chol(hess) fails; Hessian
                                        %   refreshed by fminunc when mod(isim,100)==0 or isim==1.

    % ---- ML computation sizes (ml_BVAR_*.m; defaults inside bvt.ml.*) ----
pr.ml.R = [NaN NaN 1000 NaN 1000 NaN 5000 10000];
                                        % importance-sampling draws for the integrated likelihood,
                                        %   by model: ml_BVAR_CSV.m line 10, ml_BVAR_t_CSV.m lines 13-14,
                                        %   ml_BVAR_CSV_MA.m lines 10-11, ml_BVAR_CSV_t_MA.m lines 20-21.
                                        %   Models 1/2/4/6 have analytic likelihoods (no IS).
pr.ml.nsims2 = 1000;                    % reduced-run length for the rho/psi ordinates of models 7 and 8
                                        %   (ml_BVAR_CSV_MA.m line 51, ml_BVAR_CSV_t_MA.m line 71);
                                        %   models 3/5/6 reduce over nsims sweeps instead.
pr.ml.ngrid = [NaN 700 700 700 700 300 300 299];
                                        % ordinate grid sizes by model (before inserting the starred
                                        %   value): nugrid 700 in ml_BVAR_t/ml_BVAR_t_CSV/ml_BVAR_CSV_t_MA
                                        %   (m8 line 33 keeps 700 for nu; its rho/psi grids use ngrid = 299,
                                        %   line 10), rhogrid 700 in ml_BVAR_CSV/ml_BVAR_t_CSV; ngrid = 300
                                        %   in ml_BVAR_t_MA line 11 and ml_BVAR_CSV_MA line 12; psigrid
                                        %   spans (-.99,.99) in models 4/6/7 but (-.999,.999) in model 8
                                        %   (ml_BVAR_CSV_t_MA.m line 75) - beyond the +/-.99 prior
                                        %   truncation those grid points carry lpri_psi's -1e10 penalty.

    % ---- Known legacy defects in the ML scripts (2026-09-02 audit) ----
    % Reproduced bit-for-bit by the bugcompat option of the affected bvt.ml
    % functions; the corrected computation is their default. Documented in
    % full in the affected function headers and tests/variant_map.md.
pr.ml.bug_ma_psi_line17 = true;         % ml_BVAR_MA.m line 17: the llike term s2(1)/(1+psi^2) uses the
                                        %   loop variable psi left over from the estimation run (the final
                                        %   chain draw) where psi_mean is intended (lines 11/16 use psi_mean).
                                        %   -> bvt.ml.kron_bvar_ma('bugcompat',true), fed by run_all's
                                        %   out.state.psi.
pr.ml.bug_csv_t_ma_hpsi_frozen = true;  % ml_BVAR_CSV_t_MA.m lines 42-44: the (A,Sig) ordinate loop never
                                        %   rebuilds Hpsi and never reads the stored psi draws - it consumes
                                        %   the estimation run's leftover Hpsi/psi for all nsims terms
                                        %   (ml_BVAR_CSV_MA.m lines 26-30 rebuild per draw - the intended
                                        %   pattern). -> bvt.ml.kron_bvar_csv_t_ma('bugcompat',true).
pr.ml.bug_csv_t_ma_sig_line108 = true;  % ml_BVAR_CSV_t_MA.m line 108: the reduced-run psi target uses the
                                        %   leftover last DRAW `Sig` where Sig_mean is intended (the same
                                        %   line in ml_BVAR_CSV_MA.m line 83 uses Sig_mean; this reduced
                                        %   run's own h/lam steps condition on Sig_mean via line 78).
                                        %   -> same bugcompat flag, fed by out.state.Sig.
pr.ml.quirk_csv_t_ma_intlike_first_obs = true;
                                        % intlike_BVAR_CSV_t_MA.m line 56 + its deny_h: the integrated
                                        %   likelihood treats the first TRANSFORMED observation as t with
                                        %   scale exp(h_1)*Sig - no (1+psi^2) scale factor and no
                                        %   -n/2*log(1+psi^2) constant - although the Gibbs sampler gives
                                        %   it variance (1+psi^2)*exp(h_1)*lam_1*Sig (the Gaussian
                                        %   intlike_BVAR_CSV_MA.m applies both). A likelihood-formula
                                        %   property, NOT an evaluation-point inconsistency: kept verbatim
                                        %   in BOTH bvt.ml.intlike_csv_t_ma modes and documented there.
pr.ml.quirk_csv_t_ma_lam_first_obs = true;
                                        % BVAR_CSV_t_MA.m lam step (lines 71-75) and its ml reduced run
                                        %   (ml_BVAR_CSV_t_MA.m lines 90-91) apply NO (1+psi^2) correction
                                        %   to the first observation's s2 in the lam_1 draw, while the h
                                        %   step does (line 79) and model 6 corrects its lam step
                                        %   (BVAR_t_MA.m line 97). Estimation and reduced run are mutually
                                        %   consistent; kept verbatim in both modes.

    % ---- llike_CSV_MA path hazard (2026-09-02 audit; never-merge entry in
    %      tests/variant_map.md) ----
pr.llike_csv_ma_root_has_sumh = true;   % the package-ROOT legacy/llike_CSV_MA.m line 13 includes the
                                        %   -n/2*sum(h) term; legacy/realtime_forecasts/llike_CSV_MA.m
                                        %   line 6 OMITS it. Constant in psi (cancels inside the psi-MH
                                        %   ratios and normalized grid densities at fixed h - only the
                                        %   fminunc/fminbnd floating-point path differs) but NOT as a
                                        %   likelihood ordinate. main_forecasting.m line 20 runs
                                        %   addpath('./realtime_forecasts'), after which an unqualified
                                        %   llike_CSV_MA call resolves to the REDUCED copy whenever the
                                        %   package root is not the current folder - the estimation/ML
                                        %   pipeline is protected only by cwd precedence in the legacy
                                        %   workflow. bvt.ml.llike_csv_ma pins the root semantics and all
                                        %   toolkit callers use the qualified name. (The realtime folder
                                        %   also shadows sample_h - executably identical - and sample_nu -
                                        %   log-form vs the root's normpdf-form MH ratio, mathematically
                                        %   identical, floating-point different; see bvt.sv.nu_studentt.)
end
