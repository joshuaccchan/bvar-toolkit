% chan_koop_yu2024_jbes_oisv/preset - every constant the legacy OISV FULL-SAMPLE
% ESTIMATION pipeline hard-codes (main_SVAR_fullsample.m -> func_main_SVAR_v2.m
% dispatching to the workspace scripts SVARSV_MH.m (OI) / CS_MH.m (CS)), one
% field per constant, each cited to its legacy source line. Consumed by
% run_all.m in this folder. Compiled 2026-09-02 (step 7, OISV family pass); the
% legacy folder is never modified - this file only transcribes it.
%
% The FORECAST pipeline (submain_forecasting_{OI1,OI2,CS1,CS2}.m +
% forecast_SVARSV_MH.m / forecast_CS_MH.m) is NOT functionized in this step:
% the four submain drivers are CLUSTER FRAGMENTS - their vintage loop
% `for t = T0:T-1` is commented out (submain line 52) and `t` must be assigned
% externally per cluster job (README.txt), so they are not runnable as shipped.
% Their goldens are the shipped results_mat/forecasting{OI1,OI2,CS1,CS2}-
% cluster.mat. The fragments' known numerical divergences from the estimation
% scripts (2026-09-02 audit, re-verified line-by-line from source) are recorded
% under pr.forecast as documentation-only fields so a future forecast
% functionization does not silently inherit the estimation behavior.
%
% See:
% Chan, J.C.C., Koop, G. and Yu, X. (2024). Large Order-Invariant Bayesian
% VARs with Stochastic Volatility, Journal of Business and Economic
% Statistics, 42(2): 825-837.

function pr = preset()

    % ---- data & model dimensions (func_main_SVAR_v2.m) ----
pr.data_file = 'FRED_MD_20vars.csv';        % func_main_SVAR_v2.m line 14; lives in legacy/, loaded read-only by run_all (monthly, 730 x 20)
pr.var_id1   = [4,6,12,13];                 % line 15: IP, unemployment, PCE inflation, FEDFUNDS first (the legacy "n = 23" comment there is stale - n = 20)
pr.var_id2   = [1:3,5,7:11,14:20];          % line 16: the remaining 16 series
pr.var_id    = [pr.var_id1, pr.var_id2];    % line 17; rev_option==1 (run_all flip=1) reverses it, line 20
pr.n0        = 24;                          % line 22: first 24 obs are the initial conditions Y0
pr.p         = 13;                          % line 8 (legacy comment: if p > 8, Y0/Y must change)
pr.nsim_default   = 30000;                  % main_SVAR_fullsample.m line 2
pr.burnin_default = 5000;                   % main_SVAR_fullsample.m line 3

n = numel(pr.var_id);                       % = 20; legacy gets it as [T,n] = size(Y), func_main_SVAR_v2.m line 25

    % ---- constants shared by the two estimation samplers ----
pr.ls_ridge     = .01;                      % chain-init least-squares ridge (X'X + .01*speye(k)): SVARSV_MH.m line 29, CS_MH.m line 32
pr.sv_offset    = .0001;                    % ystar = log(E.^2 + .0001): SVARSV_MH.m line 93, CS_MH.m line 93
pr.phi_init_bnd = .99;                      % chain-init truncation phi = min(phi0 + sqrt(Vphi).*randn, .99): SVARSV_MH.m line 36, CS_MH.m line 31
pr.progress_every = 2000;                   % loop-counter disp cadence: SVARSV_MH.m line 131, CS_MH.m line 131
    % horseshoe chain init (SVARSV_MH.m lines 38-42, CS_MH.m lines 40-44,
    % identical expressions in both):
    %   z_psi1 ~ 1./gamrnd(.5,1,n*p,1);  z_psi2 ~ 1./gamrnd(.5,1,(n-1)*n*p,1);
    %   z_kappa ~ 1./gamrnd(.5,1,2,1);
    %   psi_kappa1 ~ 1./gamrnd(.5,z_psi1);  psi_kappa2 ~ 1./gamrnd(.5,z_psi2);
    % run_all keeps those expressions verbatim; they carry no free constants
    % beyond the .5/1 shapes visible above.

    % ---- SVARSV_MH.m (model 1 = 'OI': order-invariant SVAR-SV) ----
pr.oi.kappa_init = [.1,.1,NaN,100];         % SVARSV_MH.m line 28: [own lag, other lag, UNUSED (NaN), intercepts]
pr.oi.nuh  = 3*ones(n,1);                   % line 9: IG shape for the SV state variances sig2
pr.oi.Sh   = .05*ones(n,1).*(pr.oi.nuh-1);  % line 9: IG scale (= .1); also drawn through at chain init, line 35
pr.oi.phi0 = .95*ones(n,1);                 % line 10: prior mean of the SV AR(1) coefficients phi
pr.oi.Vphi = .05^2*ones(n,1);               % line 10: prior variance of phi
pr.oi.B0_prior_mean = eye(n);               % line 11: Hyper.B0 - prior mean of the impact matrix rows (B0 rows ~ N(e_i, I))
pr.oi.VB0  = 1*ones(n);                     % line 11: Hyper.VB0 - prior variances of the B0 elements (all 1)
pr.oi.phi_mh_bnd = .99;                     % utility/sample_SV0para.m line 17: phi MH candidate bound |phic| < .99
                                            %   (= the bvt.sv.sv0_params DEFAULT - run_all passes nothing)
pr.oi.h_mean_in_sv = 0;                     % SVARSV_MH.m line 94: sample_SV(ystar,h,0,phi,sig2) - zero-mean log-volatilities

    % ---- CS_MH.m (model 2 = 'CS': Cholesky / triangularized SVAR-SV) ----
pr.cs.kappa_init = [.1,.1,1,100];           % CS_MH.m line 9: [own lag, other lag, alpha (declared, unused by getVbeta), intercepts]
pr.cs.beta0 = zeros(n^2*pr.p + n,1);        % line 11: Hyper.beta0 = zeros(k_beta,1), k_beta = n^2*p+n = 5220
pr.cs.alp0  = zeros(n*(n-1)/2,1);           % line 13: Hyper.alp0 - declared but never read by the sampler (the alp draw uses a zero prior mean implicitly)
pr.cs.Valp  = 1*ones(n*(n-1)/2,1);          % line 13: Hyper.Valp - prior variances of the free impact elements (all 1)
pr.cs.mu0   = zeros(n,1);                   % line 14: prior mean of the log-volatility means mu
pr.cs.Vmu   = 100*ones(n,1);                % line 14: prior variance of mu
pr.cs.nuh   = 3*ones(n,1);                  % line 15: IG shape for sig2
pr.cs.Sh    = .05*(pr.cs.nuh-1);            % line 15: IG scale (= .1; same value as pr.oi.Sh, written without the ones(n,1) factor)
pr.cs.phi0  = .95*ones(n,1);                % line 16
pr.cs.Vphi  = .05^2*ones(n,1);              % line 16
pr.cs.phi_mh_bnd = .999;                    % utility/sample_SVpara.m line 20: phi MH candidate bound |phic| < .999
                                            %   (= the bvt.sv.sv_params DEFAULT - run_all passes nothing).
                                            %   VERIFIED divergence trio: .999 (OISV sample_SVpara, CS estimation AND
                                            %   both forecast fragments) vs .99 (OISV sample_SV0para, OI) vs .998
                                            %   (chan2023_joe_mlvarsv sample_SVpara - different package, never-merge,
                                            %   see tests/variant_map.md).

    % ---- estimation-vs-forecast divergences (2026-09-02 audit) ----
    % DOCUMENTATION ONLY: the forecast fragments are not functionized in this
    % step; run_all reads nothing below this line.
pr.forecast.nsim_default   = 10000;         % submain_forecasting_*.m line 9
pr.forecast.burnin_default = 500;           % submain_forecasting_*.m line 10 (vs full-sample 5000)
pr.forecast.T0 = 108;                       % submain_forecasting_*.m line 23: first forecast origin index (legacy comment: 1970:02)
pr.forecast.vintage_loop_commented = true;  % submain line 52: `%for t = T0:T-1` - t is a CLUSTER INPUT (README.txt);
                                            %   the fragments are not runnable as shipped. Goldens = the shipped
                                            %   results_mat/forecasting{OI1,OI2,CS1,CS2}-cluster.mat.
pr.forecast.insim = 51;                     % forecast_SVARSV_MH.m line 129 / forecast_CS_MH.m line 125: inner predictive-sim draws per posterior draw
pr.forecast.horizons = 12;                  % 1-12-step-ahead, both fragments (tmpyhat 3rd dim)
pr.forecast.store_h_last_only = true;       % forecast_* store only h(end,:) (store_h_T) - no full h path stored
pr.forecast.cs_h_mu_zero = true;            % forecast_CS_MH.m line 85 passes mu = 0 to sample_SV where CS_MH.m line 94
                                            %   passes mu(ii); sample_SVpara still updates mu (line 89) and the forecast
                                            %   recursion h_Tp1 = muh + phih.*(h_Tp1-muh) + ... (lines 140/159) uses it -
                                            %   the h PATH is drawn zero-mean while mu is drawn and used downstream.
pr.forecast.oi_alpha_rewritten = true;      % forecast_SVARSV_MH.m lines 69-81 REWRITE the alpha step:
                                            %   zi = reshape(B0*(Yt-[XA(:,1:ii-1),sparse(Tt,1),XA(:,ii+1:end)])',Tt*n,1),
                                            %   Wi = kron(Xt,B0(:,ii)) (time-interleaved stacking) with an explicit
                                            %   exp(-reshape(h',Tt*n,1)) weighting - vs estimation SVARSV_MH.m lines
                                            %   76-87: yi = vec((Y-X*A)*B0')./Lambda, Wi = kron(B0(:,ii),X)./Lambda
                                            %   (equation-stacked, row-scaled). Same conditional posterior, different
                                            %   floating-point path - NOT canonicalized by bvt.samplers.eq_svar_oi
                                            %   (never-merge, tests/variant_map.md).
pr.forecast.cs_kappa_from_driver = true;    % forecast_CS_MH.m never initializes kappa (first read at line 44);
                                            %   submain_forecasting_*.m line 44 supplies kappa = [.1,.1,1,100]
                                            %   (= pr.cs.kappa_init) for model 2 - the OI branch of the submain sets
                                            %   no kappa (forecast_SVARSV_MH.m line 20 sets its own [.1,.1,NaN,100]).
pr.forecast.phi_bnds_unchanged = true;      % both fragments call the same utility/sample_SVpara (.999) and
                                            %   utility/sample_SV0para (.99) as estimation - no forecast-side divergence.
pr.forecast.no_seed_lines = true;           % neither fragment carries a clock-seed line (forecast_CS_MH has none at all;
                                            %   CS_MH.m line 49 carries the estimation line COMMENTED OUT; the only ACTIVE
                                            %   clock-seed line in the package is SVARSV_MH.m line 24).
end
