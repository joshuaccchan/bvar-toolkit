% chan2023_jbes_hybtvp/preset - every constant of main_HYB_TVPSV.m in one place,
% each tagged with the legacy line it comes from. run_all reads this; nothing
% here is computed from the data.
%
% Chan, J.C.C. (2023). Large Hybrid Time-Varying Parameter VARs, Journal of
% Business and Economic Statistics, 41(3): 890-905.

function pr = preset()

    % ---- run length and data ----
pr.nsim    = 50000;                         % main_HYB_TVPSV.m 16
pr.burnin  = 1000;                          % 17
pr.p       = 2;                             % 20
pr.data_file  = 'macrodata_Q_2018Q4.csv';   % 21
pr.data_range = 'A1:IN238';                 % 21: 1959Q1-2018Q4
pr.progress_every = 1000;                   % 181: loop-counter disp cadence

    % ---- variable selections (line 23 is the active one; 22 and 24 are
    % commented out in the legacy file and kept here for the other dimensions) ----
pr.varid_n3  = [1,95,59];                                       % 22
pr.varid     = [1,95,59,144,22,133];                            % 23, ACTIVE (n = 6)
pr.varid_n20 = [1,95,59,144,22,133,160,2,18,23,35,57,76,97,120,123,135,138,148,161];  % 24

    % ---- switches ----
pr.is_gamfixed   = false;                   % 33
pr.is_kappafixed = false;                   % 34

    % ---- prior ----
pr.kappa0 = [.4 .04^2];                     % 37: initial (kappa1, kappa2)
    % bvar.priors.vtheta takes kappa as a 4-vector; this package hard-codes the
    % last two inside its own getVtheta (kappa_3 = .2 for the impact matrix,
    % kappa_4 = 1 for the intercepts). See that function's header.
pr.kappa_3 = .2;
pr.kappa_4 = 1;
pr.nuh0_scalar = 3;                         % 41: nuh0 = 3*ones(n,1)
pr.Sh0_factor  = .1;                        % 41: Sh0 = .1*(nuh0-1)
pr.Vsigbeta_offdiag = .01^2;                % 42
pr.Vsigbeta_intercept = .1^2;               % 43
pr.Vsigalp = .01^2;                         % 44
pr.c01 = [1, 1/.04];                        % 45: prior for kappa1
pr.c02 = [1, 1/.04^2];                      % 46: prior for kappa2
pr.ah_scalar = 0;                           % 47: ah = zeros(n,1)
pr.Vh_scalar = 10;                          % 47: Vh = 10*ones(n,1)
pr.ap = .5*ones(1,2);                       % 48
pr.bp = .5*ones(1,2);                       % 48
pr.p0_init = .5;                            % initialize.m 25: fixed at 0.5 so the
                                            % gam prior is symmetric and does not
                                            % affect the ML value
pr.sv_offset_init = 1e-4;                   % initialize.m 35: log(U^2 + 1e-4) at
                                            % initialization ONLY - the main loop
                                            % (main 118) uses log(Ui.^2) with no
                                            % offset

    % ---- reporting ----
pr.eqtext = ["real GDP", "PCE inflation", "Unemployment", "Fed funds rate",...
    "Industrial production index", "Real average hourly earnings in manufacturing"...
    "M1", "Real PCE", "Real disposable personal income", "Industrial production: final products", ...
    "All employees: total nonfarm", "Civilian employment", ...
    "Nonfarm business section: hours of all persons", "GDP deflator",...
    "CPI", "PPI", "Nonfarm business sector: real compensation per hour",...
    "Nonfarm business section: real output per hour",...
    "10-year treasury constant maturity rate", "M2"];   % 203-210, verbatim. The
                                            % comma after "...in manufacturing" is
                                            % missing in the legacy file; inside [ ]
                                            % a space separates too, so this is
                                            % still 20 entries and the labels line
                                            % up (checked, not assumed).

    % ---- core functions this package's driver uses ----
pr.core_used = {'bvar.priors.resid_var_allvars_ridge (get_resid_var_v2)', ...
    'bvar.priors.minnesota_C (get_C)', 'bvar.priors.vtheta (getVtheta, with kappa(3:4) supplied)', ...
    'bvar.samplers.eq_hyb_tvp (sample_gam_thetai_ver2)', ...
    'bvar.sv.ksc_rw_h0 (sample_SVRW; bitwise-verified)', ...
    'bvar.util.gam_mode (get_gammode)', 'third_party/gigrnd.m'};

    % ---- the legacy clock seeding, recorded for the equivalence test ----
pr.notes.clock_seed_line = 'main_HYB_TVPSV.m:85';   % randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
pr.notes.dead_vtheta_call = true;           % main 40 computes [Valp,Vbeta] before
                                            % the loop; line 91 recomputes them as
                                            % the first act of every sweep, so the
                                            % line-40 values are never read
end
