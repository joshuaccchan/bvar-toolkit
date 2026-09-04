% chan2023_joe_mlvarsv/preset - every constant the legacy ml_varsv estimation
% pipeline hard-codes (main_varsv.m dispatching to the workspace scripts
% VAR_NCP.m, VAR_CSV.m, VAR_ARSV_redu.m, VAR_FSV.m, VAR_ARSVO_redu.m), one
% field per constant, each cited to its legacy source line. Consumed by
% run_all.m in this folder. Compiled 2026-09-03 (step 9); the legacy folder is
% never modified - this file only transcribes it.
%
% The marginal-likelihood routines (utility/ml_var_csv.m, ml_var_arsv_redu.m,
% ml_var_fsv.m, ml_var_arsvo_redu.m) are a separate phase and are not
% functionized here; the constants they own (M, flag_marg) are recorded under
% pr.ml as documentation only. Model 1 is the exception: VAR_NCP.m computes its
% log marginal likelihood inline and analytically, so run_all reproduces it.
%
% Hyperparameters whose length depends on n are stored as the scalar the legacy
% writes, with the expansion rule in the comment; run_all expands them once n is
% known (n = numel(varid)).
%
% See:
% Chan, J.C.C. (2023). Comparing stochastic volatility specifications for large
% Bayesian VARs, Journal of Econometrics, 235(2), 1419-1446.

function pr = preset()

    % ---- data, design, run length (main_varsv.m) ----
pr.data_file = 'macrodata_Q_2019Q4.csv';    % line 33; lives in legacy/, read-only (242 x 248, quarterly)
pr.varid     = [1,22,59,120,133,144,148,2,35,57,81,95,152,160,245];   % line 35: the active n = 15 selection
pr.varid_n7  = [1,22,59,120,133,144,148];   % line 34, commented out: the n = 7 selection
pr.varid_n30 = [1,22,59,120,133,144,148,2,35,57,81,95,152,160,245, ...
    3,18,23,37,58,76,83,97,123,138,145,147,157,161,199];              % lines 36-37, commented out: n = 30
pr.n0 = 8;                                  % line 40: Y0 = data(1:8,:) (the "first 9 obs" comment there is stale)
pr.p  = 4;                                  % line 24 (legacy comment: if p > 8, Y0/Y must change)
pr.r  = 2;                                  % line 25: number of latent factors, VAR-FSV only
pr.nsim_default   = 10000;                  % line 26
pr.burnin_default = 1000;                   % line 27
pr.progress_every = 10000;                  % loop-counter disp cadence: VAR_CSV.m 81, VAR_ARSV_redu.m 104, VAR_FSV.m 110, VAR_ARSVO_redu.m 126

    % ---- Minnesota / natural-conjugate shrinkage constants (main_varsv.m) ----
pr.kappa  = .2^2;                           % line 56: VAR-NCP / VAR-CSV single shrinkage (VAR-CSV then samples it)
pr.kappa1 = .2^2;                           % line 57: own-lag shrinkage (models 3/4/5)
pr.kappa3 = 100;                            % line 57: intercept prior variance (all models)
pr.kappa4 = .2^2;                           % line 57: impact-matrix (B0) shrinkage, models 3/5
pr.kappa2 = (.2^2)^2;                       % line 61: cross-lag shrinkage when is_kappasym is false
                                            %   line 59: is_kappasym true -> kappa2 = kappa1 (= pr.kappa1)

    % ---- constants shared by the SV samplers ----
pr.sv_offset  = .0001;                      % ystar = log(u^2 + .0001): VAR_ARSV_redu.m 79, VAR_ARSVO_redu.m 86;
                                            %   VAR_FSV.m 76 writes the same value as 1e-4
pr.phi_mh_bnd = .998;                       % utility/sample_SVpara.m line 22: phi MH candidate bound |phic| < .998.
                                            %   Passed explicitly to bvar.sv.sv_params, whose default is the OISV .999.
pr.nuh_scalar = 3;                          % IG shape for sig2: main_varsv.m 76 (scalar), 89/118 (3*ones(n,1)), 102 (3*ones(n+r,1))
pr.Sh_factor  = .1;                         % IG scale = Sh_factor*(nuh-1) = .2: main_varsv.m 76/89/102/118
pr.mu0_scalar = 0;                          % prior mean of the log-vol means: main_varsv.m 90/119 (zeros(n,1)), 103 (zeros(n+r,1))
pr.Vmu_scalar = 100;                        % prior variance of mu: main_varsv.m 90/119/103
pr.phi0_scalar = .98;                       % prior mean of phi: main_varsv.m 77 (scalar), 91/120 (n), 104 (n+r)
pr.Vphi_scalar = .05^2;                     % prior variance of phi: main_varsv.m 77/91/104/120

    % ---- VAR-NCP (model 1, VAR_NCP.m: analytic, no MCMC) ----
    % Prior from prior_NCP(p,kappa,kappa3,Y0,Y) [main_varsv.m 67]
    % = bvar.priors.niw(p,[kappa kappa3],Y0,Y,'mlvarsv_ncp'); nu0 = n+3,
    % S0 = diag(sig2) are set inside that constructor (prior_NCP.m line 35).
pr.ncp.cp_ml_inline = true;                 % VAR_NCP.m lines 18-21: the log-ML is inline and analytic - in scope here

    % ---- VAR-CSV (model 2, VAR_CSV.m) ----
    % Prior from prior_NCP(p,kappa,kappa3,Y0,Y) [main_varsv.m 74], re-derived at
    % VAR_CSV.m 18 and re-derived every sweep at line 36 with the current kappa.
pr.csv.c0 = [1,1/.2^2];                     % main_varsv.m 75: GIG hyperparameters for the kappa draw
pr.csv.C_from_VA = true;                    % VAR_CSV.m 19: C = Hyper.VA/kappa - the kappa-free part of VA
pr.csv.h_init_forced_accept = true;         % VAR_CSV.m 27: sample_CSV(s2,phi,sig2,zeros(T,1),n,true) at chain init
pr.csv.forced_accept_before = 20;           % VAR_CSV.m 52: `if isim < 20` forces the h proposal through during early burn-in
pr.csv.kappa_shape_offset = 'n^2*p/2';      % VAR_CSV.m 69: gigrnd(c0(1)-n^2*p/2, 2*c0(2), tmpc, 1) over the k-1 non-intercept rows
pr.csv.plot_T_id = [1961, 2019.75];         % VAR_CSV.m 102: linspace for the exp(h/2) figure - cosmetic, not reproduced by run_all

    % ---- VAR-SV (model 3, VAR_ARSV_redu.m) ----
    % Priors: prior_Minn(p,kappa1,kappa2,kappa3,Y0,Y) [main_varsv.m 86]
    % = bvar.priors.minn(...,4); prior_B0(Y0,Y,kappa4) [87] = bvar.priors.impact_B0.
    % Both are re-derived every sweep (lines 42 and 61) with the current kappas.
pr.sv.c0 = [1,1/.2^2; 1,1/.2^2; 1,1];       % main_varsv.m 88: GIG rows for kappa1, kappa2, kappa4
pr.sv.C_beta_from_Vbeta = true;             % VAR_ARSV_redu.m 35: C_beta = Hyper.Vbeta/kappa4
pr.sv.kappa4_always_drawn = true;           % VAR_ARSV_redu.m 88-89 sit before the is_kappafixed branch: kappa4 is
                                            %   drawn every sweep even when kappa1/kappa2 are held fixed
pr.sv.kappa4_shape_offset = 'n*(n-1)/4';    % line 89: gigrnd(c0(3,1)-n*(n-1)/4, 2*c0(3,2), tmpc4, 1)
pr.sv.kappa1_shape_offset = 'n*p/2';        % line 100 (asymmetric branch)
pr.sv.kappa2_shape_offset = '(n-1)*n*p/2';  % line 101
pr.sv.kappa_sym_shape_offset = 'n^2*p/2';   % line 95: the is_kappasym branch pools both index sets into one draw

    % ---- VAR-FSV (model 4, VAR_FSV.m) ----
    % Prior: prior_Minn only [main_varsv.m 100]; no impact matrix, no kappa4.
pr.fsv.c0 = [1,1/.2^2; 1,1/.2^2;];          % main_varsv.m 101: GIG rows for kappa1, kappa2
pr.fsv.l0 = 0;                              % main_varsv.m 105: prior mean of the free factor loadings
pr.fsv.Vl = 1;                              % main_varsv.m 106: prior variance of the free factor loadings
pr.fsv.A_init_VA = [100, .1];               % VAR_FSV.m 19: ridge 1./[100, .1*ones(1,k-1)] on the chain-init least squares
pr.fsv.L_init_lower = 1;                    % VAR_FSV.m 20: L = [eye(r); ones(n-r,r)]
pr.fsv.phi_init = [.9, .09];                % VAR_FSV.m 24: phi = .9 + .09*rand(n+r,1)
pr.fsv.sig2_init = [.05, .05];              % VAR_FSV.m 25: sig2 = .05 + .05*rand(n+r,1)
pr.fsv.mu_init_halved = true;               % VAR_FSV.m 23: mu = [log(var(Y)'/2); log(mean(var(Y)'))*ones(r,1)]
pr.fsv.h_init_factor_rows = true;           % VAR_FSV.m 27: h = [zeros(T,n) repmat(mu(n+1:end)',T,1)], then the first
                                            %   n columns are overwritten by getARh_approx1N (line 29)
pr.fsv.store_A_is_running_sum = true;       % VAR_FSV.m 10/103/118: store_A accumulates A and is divided by nsim - no per-draw store

    % ---- VAR-SVO (model 5, VAR_ARSVO_redu.m) ----
    % Priors as VAR-SV [main_varsv.m 115-120] plus the outlier block.
pr.svo.c0 = [1,1/.2^2; 1,1/.2^2; 1,1];      % main_varsv.m 117 (identical to pr.sv.c0)
pr.svo.p0a = 10/4;                          % main_varsv.m 121: beta prior shape a for the outlier probability po
pr.svo.p0b = (1-1/16)*40;                   % main_varsv.m 121: beta prior shape b
pr.svo.ngrid = 31;                          % VAR_ARSVO_redu.m 8: number of outlier-size grid points above 1
pr.svo.o_grid_range = [2, 20];              % VAR_ARSVO_redu.m 9: o_grid = [1; linspace(2,20,ngrid)'] (32 points total)
pr.svo.po_init = 1/16;                      % VAR_ARSVO_redu.m 22
pr.svo.o_init = 1;                          % VAR_ARSVO_redu.m 23: o = ones(T,1)

    % ---- marginal likelihood (run_ml.m in this folder; step 10) ----
    % run_all reads only pr.ml.cp_ml (model 1); run_ml reads M and flag_marg.
pr.ml.M = 10000;                            % main_varsv.m 28: importance-sampling draws
pr.ml.cp_ml = 1;                            % main_varsv.m 29
pr.ml.flag_marg = 2;                        % main_varsv.m 30: integrate out A/(A,Sig) and sig2.
                                            %   Only ml_var_fsv.m implements flag_marg = 1; the other three define
                                            %   prior/gIS in the case-2 branch alone
pr.ml.routines = {'(inline in VAR_NCP.m)', 'ml_var_csv', 'ml_var_arsv_redu', ...
    'ml_var_fsv', 'ml_var_arsvo_redu'};     % main_varsv.m dispatch -> bvar.ml.mlvarsv_* (step 10)
pr.ml.getISden_ARSS_used_by_ml_only = true; % utility/getISden_ARSS.m has no caller in the five estimation scripts;
                                            %   it is read only by the ml_var_* routines (glob-verified 2026-09-03)
pr.ml.kappa3 = 100;                         % re-hard-coded INSIDE every ml_var_* routine (ml_var_csv.m 9,
                                            %   ml_var_arsv_redu.m 14, ml_var_arsvo_redu.m 14, ml_var_fsv.m 17),
                                            %   duplicating pr.kappa3; kept there, not passed in
pr.ml.nbatch = 50;                          % M is rounded up to a multiple of 50 and the weights are averaged in
                                            %   50 batches; the reported standard error is std(batch)/sqrt(50)
pr.ml.h_is_rho_bnd = [-.99, .99];           % getISden_ARSS.m 12: fminbnd range for the IS density's AR coefficient
pr.ml.phi_is_trunc = [-1, 1];               % every routine draws and scores phi on (-1,1), while the SAMPLER
                                            %   truncates its MH candidate at |phi| < .998 (pr.phi_mh_bnd). The IS
                                            %   estimator stays valid - the (-1,1) prior is what it integrates - but
                                            %   the posterior the draws come from is the .998 one. Family-wide, unchanged
pr.ml.svo_defects = {'o prior mass over 32 atoms, not the sampler''s 31 (ml_var_arsvo_redu.m 180)', ...
    'o_hat(o_idx) linear index -> column 1 for T >= 32 (line 181)', ...
    'c1 omits -n*sum(log(o)), the outlier Jacobian (line 153)'};
                                            % bvar.ml.mlvarsv_arsvo_redu corrects all three by default and reproduces
                                            %   the legacy under 'bugcompat', true. tests/variant_map.md has the audit
                                            %   and the effect on the published value
pr.ml.fsv_2pi_cancels = true;               % ml_var_fsv.m 125/159: lh_prior and lh_g both drop T(n+r)/2*log(2*pi);
                                            %   the omissions cancel in llike + lh_prior - lh_g. Not a defect
pr.ml.fsv_big_sig2_unused_at_flag2 = true;  % ml_var_fsv.m 67-70 draws big_sig2 even under flag_marg = 2, where
                                            %   nothing reads it - rng-consuming, numerically inert
pr.ml.arsv_dead_gamfit = true;              % ml_var_arsv_redu.m 47-51 / ml_var_arsvo_redu.m 47-51 fit a gamma to
                                            %   1./sig2 and never use it (the lines that would are commented out)

    % ---- divergences a reader should not trip over ----
pr.notes.clock_seed_lines = {'VAR_CSV.m:32', 'VAR_ARSV_redu.m:38', ...
    'VAR_FSV.m:34', 'VAR_ARSVO_redu.m:45'};  % all four active: randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
                                             %   run_all does not reproduce them (they are irreproducible and switch
                                             %   MATLAB to the v4/v5 generators); VAR_NCP.m has none.
pr.notes.is_kappasym_ignored_by_csv = true;  % main_varsv.m 68-73: model 2's name/branch ignores is_kappasym; only
                                             %   is_kappafixed gates its kappa draw (VAR_CSV.m 66)
pr.notes.alp_row_vs_column = true;           % bvar.samplers.alp_tri_cs returns a 1 x k_beta row; VAR_ARSV_redu.m 31
                                             %   keeps beta as a k_beta x 1 column and both store_beta(isave,:) = beta'
                                             %   (line 111) and tmpc4 = sum(beta.^2./C_beta) (line 88) depend on that
                                             %   orientation - run_all transposes at the call site
pr.notes.dead_U_prealloc = true;             % VAR_ARSV_redu.m 43 / VAR_ARSVO_redu.m 50: U = zeros(T,n) is never read
                                             %   before being overwritten; kept in run_all as a dead assignment
pr.notes.nuub_absent = true;                 % there is no constant named nuub (or nu_ub / nub) anywhere in this
                                             %   legacy folder - grep-verified 2026-09-03. A carried-over lead that
                                             %   does not correspond to this package.
end
