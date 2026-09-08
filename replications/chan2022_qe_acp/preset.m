% chan2022_qe_acp/preset - every constant of main_ACP_apps.m and
% main_ACP_jointden.m in one place, each tagged with the legacy line it comes
% from. run_all and run_jointden read this; nothing here is computed from data.
%
% Chan, J.C.C. (2022). Asymmetric Conjugate Priors for Large Bayesian VARs,
% Quantitative Economics, 13(3): 1145-1169.

function pr = preset()

    % ---- shared ----
pr.p = 5;                                   % main_ACP_apps.m 16 (p > 8 needs Y0/Y changed)
pr.n_init = 8;                              % 40: first 8 observations are initial conditions
pr.data_file = 'database_2019Q4.xlsx';      % 28
pr.horizon = 36;                            % 26: steps for the impulse responses
pr.nbatch = 50000;                          % 25: posterior draws sampled per batch
pr.progress_every = 1e6;                    % 140: report cadence, in total draws

    % ---- dataset 1: 6 variables, 5 shocks ----
pr.d1.var_id = 1:6;                         % 32: GDP, deflator, interest rate, investment, S&P, spread
pr.d1.idx_ns = [1,2,4,5];                   % 33: variables entering in levels
pr.d1.nsim = 5000;                          % 19: accepted draws to collect
pr.d1.kappa = [1, 1, 1, 100];               % 50: fixed, not optimized
    % 62-67, built the legacy way: one named column vector per shock, in the
    % order supply, demand, monetary, investment, financial
pr.d1.S = [[1,-1,NaN,NaN,1,NaN]', ...       % supply
           [1,1,1,NaN,NaN,NaN]', ...        % demand
           [1,1,-1,NaN,NaN,NaN]', ...       % monetary
           [1,1,1,NaN,-1,NaN]', ...         % investment
           [1,1,1,NaN,1,NaN]'];             % financial
pr.d1.Rineq = [-1,0,0,1,0,0; 1,0,0,-1,0,0; 1,0,0,-1,0,0];   % 70
pr.d1.Ridx = [2,4,5];                       % 71: columns Rineq applies to
pr.d1.titles = ["GDP" "GDP Deflator" "3-month Tbill" "Investment" "S&P 500" "Spread"];  % 153

    % ---- dataset 2: 15 variables, 5 shocks ----
pr.d2.var_id = 1:15;                        % 37
pr.d2.idx_ns = [1,2,4,5,10,11,12,13,15];    % 38
pr.d2.nsim = 1000;                          % 21: "might take a few days"
pr.d2.kappa_init = [.04,.0016];             % 53: starting values for get_OptKappa
    % 74-79, same construction at n = 15
pr.d2.S = [[1,-1,NaN,NaN,1,NaN,NaN,NaN,NaN,-1,-1,NaN,1,NaN,1]', ...    % supply
           [1,1,1,NaN,NaN,NaN,NaN,NaN,NaN,1,1,NaN,1,1,NaN]', ...       % demand
           [1,1,-1,NaN,NaN,NaN,NaN,NaN,NaN,1,1,NaN,1,-1,NaN]', ...     % monetary
           [1,1,1,NaN,-1,NaN,NaN,NaN,NaN,1,1,NaN,1,1,-1]', ...         % investment
           [1,1,1,NaN,1,NaN,NaN,NaN,NaN,1,1,NaN,1,1,1]'];              % financial
pr.d2.Rineq = [-1,0,0,1,zeros(1,11); 1,0,0,-1,zeros(1,11); 1,0,0,-1,zeros(1,11)];  % 82
pr.d2.Ridx = [2,4,5];                       % 83
pr.d2.varidx_prices = [2,10,11];            % 182: the second figure's variables
pr.d2.titles_prices = ["GDP Deflator" "CPI" "PCE"];                                % 181

    % ---- figure limits (main_ACP_apps.m 150-151), kept for a caller that plots ----
pr.ylim_u = [.008, .003, .4, .03, .015, .2];
pr.ylim_l = [-.001, -.001, -.2, -.01, -.005, -.4];

    % ---- main_ACP_jointden.m ----
pr.jd.idx_ns = [1,2,4,5,10,11,12,13,15];    % main_ACP_jointden.m 20
pr.jd.var_id = 1:15;                        % 21-22
pr.jd.kappa = [.04,.04^2,1,100];            % 31
pr.jd.kappa1_grid = 0.01:.001:.2;           % 34
pr.jd.kappa2_grid = .001:.0002:0.012;       % 34
pr.jd.subjective = [.04,.0016];             % 58: the point marked on the contour

    % ---- notes ----
pr.notes.xlsread_required = true;           % readmatrix returns a DIFFERENT SHAPE on
                                            % this file ([141 16] against xlsread's
                                            % [140 15]): it keeps a header row and an
                                            % index column that xlsread's numeric mode
                                            % drops. The legacy call is kept.
pr.notes.no_clock_seed = true;              % neither driver seeds from the clock, so
                                            % the equivalence test needs no seed patch
pr.core_used = {'bvar.priors.resid_var_ar4 (get_resid_var)', ...
    'bvar.priors.acp_redu / acp_stru (prior_ACP_redu / prior_ACP_stru)', ...
    'bvar.priors.acp_opt_kappa (get_OptKappa, and get_OptSymKappa under ''symmetric'')', ...
    'bvar.samplers.acp_theta_sig (sample_ThetaSig)', ...
    'bvar.structural.reduced_form (getReducedForm)', ...
    'bvar.structural.qr_sign (QR)', ...
    'bvar.structural.sign_restrict (the inline block, main_ACP_apps.m 103-131)', ...
    'bvar.structural.irf_redu (IRredu)', 'bvar.ml.acp (ml_VAR_ACP)'};
end
