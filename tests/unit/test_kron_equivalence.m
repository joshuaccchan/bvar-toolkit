function test_kron_equivalence
% seeded draw-for-draw equivalence of the functionized Kronecker pipeline
% (replications/chan2020_jbes_kronecker/run_all.m + run_ml.m + bvt.ml.* +
% bvt.sv.csv_armh/nu_studentt + bvt.priors.niw('kron_script') +
% bvt.util.build_lags) with the legacy workspace scripts of Chan (2020,
% JBES), for ALL EIGHT models, estimation AND marginal likelihood in one
% continuous stream (the legacy cp_ml = 1 pipeline: main_BVAR.m primes the
% workspace, the estimation script runs, its tail dispatches ml_BVAR_* in
% the same workspace), at small nsim from tempdir copies. Asserts isequal on
% all stores, acceptance counters, script-tail summaries, every ML piece
% (ML, llike, lpri, lpost, the final store_lpost), and the terminal rng
% state. Models 4 and 8 run the legacy comparison under
% run_ml(...,'bugcompat',true) - the legacy ml scripts' leftover-workspace
% defects reproduced bitwise - and then assert the CORRECTED default mode
% differs exactly where each bug lives and matches everywhere else.
%
% SOLE PATCH to the legacy scripts: each of the seven MCMC estimation
% scripts carries the clock-seed line
%   randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
% ACTIVE (asserted: exactly one occurrence each, not commented) - it
% re-seeds the global stream from the wall clock (irreproducible) AND
% switches MATLAB to the legacy v4/v5 generators, so it is removed from the
% tempdir copies. BVAR.m (model 1, no MCMC) and all seven ml_BVAR_*.m, four
% intlike_*.m and the helper files are asserted seed-line-free and copied
% BYTE-VERBATIM. Model 2/5/6/8 chain-init gamrnd draws sit before the
% removed line, so seeding rng(seed,'twister') before dispatch aligns the
% whole run, chain init included. The harness replicates main_BVAR.m lines
% 26-31 (data/dims priming) and overrides its hard-coded nsims/burnin, which
% the scripts read from the workspace.
root = getappdata(0, 'bvt_repo_root');
leg = fullfile(root, 'replications', 'chan2020_jbes_kronecker', 'legacy');
repdir = fullfile(root, 'replications', 'chan2020_jbes_kronecker');

nsims = 20; burnin = 10; seed = 20260902;
    % runtime is dominated by the ml_BVAR_CSV_MA / ml_BVAR_CSV_t_MA reduced
    % runs (hard-coded nsims2 = 1000) and intlike R (5000/10000): ~26 s and
    % ~34 s per side; the whole test runs in roughly 3 minutes.

% --- tempdir: patched MCMC scripts + byte-verbatim everything else ---
tmp = tempname; mkdir(tmp);
ctmp = onCleanup(@() cleanup_tmp(tmp));   % rmpath BEFORE rmdir, warning-free
seedline = 'randn(''seed'',sum(clock*100)); rand(''seed'',sum(clock*1000));';

mcmc_scripts = {'BVAR_t.m','BVAR_CSV.m','BVAR_MA.m','BVAR_t_CSV.m', ...
    'BVAR_t_MA.m','BVAR_CSV_MA.m','BVAR_CSV_t_MA.m'};
for kf = 1:numel(mcmc_scripts)
    txt = fileread(fullfile(leg, mcmc_scripts{kf}));
    assert(numel(strfind(txt, seedline)) == 1, ...
        'expected exactly one clock-seed line in %s', mcmc_scripts{kf});
    assert(isempty(strfind(txt, ['%' seedline])), ...
        'the clock-seed line in %s is expected to be ACTIVE', mcmc_scripts{kf});
    txt = strrep(txt, seedline, ...
        '% [clock-seed line removed by test_kron_equivalence - the sole patch]');
    fid = fopen(fullfile(tmp, mcmc_scripts{kf}), 'w');
    fwrite(fid, txt);
    fclose(fid);
end

verbatim = {'BVAR.m', ...
    'ml_BVAR_t.m','ml_BVAR_CSV.m','ml_BVAR_MA.m','ml_BVAR_t_CSV.m', ...
    'ml_BVAR_t_MA.m','ml_BVAR_CSV_MA.m','ml_BVAR_CSV_t_MA.m', ...
    'intlike_BVAR_CSV.m','intlike_BVAR_t_CSV.m','intlike_BVAR_CSV_MA.m', ...
    'intlike_BVAR_CSV_t_MA.m', ...
    'construct_prior_A.m','sample_h.m','sample_nu.m', ...
    'llike_MA.m','llike_CSV_MA.m','lniwpdf.m','linvgammpdf.m'};
for kf = 1:numel(verbatim)
    txt = fileread(fullfile(leg, verbatim{kf}));
    assert(isempty(strfind(txt, seedline)), ...
        '%s unexpectedly carries a clock-seed line', verbatim{kf});
    copyfile(fullfile(leg, verbatim{kf}), fullfile(tmp, verbatim{kf}));
end

addpath(tmp);   % removed by cleanup_tmp via ctmp, before the folder is deleted
addpath(repdir); cp2 = onCleanup(@() rmpath(repdir)); %#ok<NASGU>

% three packages define run_all - pin the resolution to this one
resolved = which('run_all');
assert(strncmpi(resolved, repdir, numel(repdir)), ...
    'run_all must resolve from the kronecker package, got %s', resolved);
% the ROOT llike_CSV_MA must win inside the legacy dispatch (path hazard)
resolved = which('llike_CSV_MA');
assert(strncmpi(resolved, tmp, numel(tmp)), ...
    'llike_CSV_MA must resolve from the tempdir ROOT copy, got %s', resolved);

C4 = []; C8 = [];
for model = 1:8
    L = run_legacy(model, leg, tmp, nsims, burnin, seed);
    C = run_ml(model, nsims, burnin, seed, 'bugcompat', true);
    sC = rng;

    switch model
        case 1
            est = {'Ahat','Shat','KA'}; ml = {'llike'};
        case 2
            est = {'store_A','store_Sig','store_nu','store_lam','store_pnu', ...
                'countnu','A_mean','Sig_mean','nu_mean','pnu_mean'};
            ml = {'llike','lpri','lpost','store_lpost'};
        case 3
            est = {'store_A','store_Sig','store_h','store_theta','counth', ...
                'A_mean','Sig_mean','h_mean','theta_mean','CSV_std_mean'};
            ml = {'llike','lpri','lpost','store_lpost'};
        case 4
            est = {'store_A','store_Sig','store_psi','store_ppsi','countpsi', ...
                'A_mean','Sig_mean','psi_mean'};
            ml = {'llike','lpri','lpost','store_lpost','den_psi'};
        case 5
            est = {'store_A','store_Sig','store_h','store_lam','store_theta', ...
                'store_pnu','counth','countrho','countnu', ...
                'A_mean','Sig_mean','theta_mean','h_mean','pnu_mean','CSV_std_mean'};
            ml = {'llike','lpri','lpost','store_lpost'};
        case 6
            est = {'store_A','store_Sig','store_lam','store_theta','store_pnu', ...
                'store_ppsi','countnu','countpsi', ...
                'A_mean','Sig_mean','theta_mean','ppsi_mean','pnu_mean'};
            ml = {'llike','lpri','lpost','store_lpost'};
        case 7
            est = {'store_A','store_Sig','store_h','store_theta','store_ppsi', ...
                'counth','countrho','countpsi', ...
                'A_mean','Sig_mean','h_mean','theta_mean','ppsi_mean','CSV_std_mean'};
            ml = {'llike','lpri','lpost','store_lpost'};
        case 8
            est = {'store_A','store_Sig','store_h','store_lam','store_theta', ...
                'store_pnu','store_ppsi','count_h','count_rho','count_psi','count_nu', ...
                'A_mean','Sig_mean','CSV_std_mean','theta_mean','ppsi_mean','pnu_mean'};
            ml = {'llike','lpri','lpost','store_lpost'};
    end
    for kf = 1:numel(est)
        assert(isequal(L.(est{kf}), C.(est{kf})), ...
            'model %d: %s differs', model, est{kf});
    end
    for kf = 1:numel(ml)
        assert(isequal(L.(ml{kf}), C.ml.(ml{kf})), ...
            'model %d: ml %s differs', model, ml{kf});
    end
    assert(isequal(L.ML, C.ML), 'model %d: ML differs', model);
    if model == 3
        % the legacy ml continues the ESTIMATION countrho counter
        assert(isequal(L.countrho, C.ml.countrho), 'model 3: final countrho differs');
    end
    assert(isequal(L.rngstate, sC.State), 'model %d: rng call sequence differs', model);

    if model == 4, C4 = C; end
    if model == 8, C8 = C; end
end

% --- corrected default mode: differs exactly where each legacy bug lives ---
% model 4 (ml_BVAR_MA.m line 17 leftover psi): only the llike term moves
C4c = run_ml(4, nsims, burnin, seed);   % default = corrected
assert(~isequal(C4c.ml.llike, C4.ml.llike), ...
    'model 4 corrected: llike should differ from bugcompat (the line-17 term)');
assert(~isequal(C4c.ML, C4.ML), 'model 4 corrected: ML should differ from bugcompat');
assert(isequal(C4c.ml.lpri, C4.ml.lpri) && isequal(C4c.ml.lpost, C4.ml.lpost) ...
    && isequal(C4c.ml.den_psi, C4.ml.den_psi) ...
    && isequal(C4c.ml.store_lpost, C4.ml.store_lpost), ...
    'model 4 corrected: only the llike piece should move');
assert(isequal(C4c.ml.psi_llike, C4c.ml.psi_mean) ...
    && isequal(C4.ml.psi_llike, C4.state.psi), ...
    'model 4: llike evaluation psi must be psi_mean (corrected) / final chain draw (bugcompat)');

% model 8 (frozen Hpsi/psi ordinate loop + line-108 leftover Sig): the
% intlike (drawn first, same stream position, same inputs) is unchanged; the
% (A,Sig) ordinate moves; the sigh2/nu ordinates are untouched; the reduced
% run (and hence ML) diverges
C8c = run_ml(8, nsims, burnin, seed);   % default = corrected
assert(isequal(C8c.ml.llike, C8.ml.llike), ...
    'model 8 corrected: the integrated likelihood should be unchanged');
assert(~isequal(C8c.ml.lpost(1), C8.ml.lpost(1)), ...
    'model 8 corrected: the (A,Sig) ordinate should differ (frozen-Hpsi bug)');
assert(isequal(C8c.ml.lpost(2:3), C8.ml.lpost(2:3)), ...
    'model 8 corrected: the sigh2/nu ordinates should be unchanged');
assert(~isequal(C8c.ML, C8.ML), 'model 8 corrected: ML should differ from bugcompat');
end

% -------------------------------------------------------------------------
function cleanup_tmp(tmp)
% deterministic teardown order: path entry first, then the folder itself
if any(strcmpi(strsplit(path, pathsep), tmp))
    rmpath(tmp);
end
if exist(tmp, 'dir')
    rmdir(tmp, 's');
end
end

% -------------------------------------------------------------------------
function L = run_legacy(model, leg, tmp, nsims, burnin, seed) %#ok<INUSD> % nsims/burnin/cp_ml are read by the dispatched scripts from this workspace
% prime the workspace exactly as main_BVAR.m lines 26-31 do, dispatch the
% model's estimation script from the tempdir (which chains into its
% ml_BVAR_* companion because cp_ml = 1), and collect everything compared
cp_ml = 1; %#ok<NASGU>
p = 4; %#ok<NASGU>
data_Q = load(fullfile(leg, 'data_Q.csv'));             % legacy data, read-only
data = data_Q(:, [1:3 6:15 17 19:24]);
Y0 = data(1:4, :); %#ok<NASGU>
shortY = data(5:end, :);
[T, n] = size(shortY); %#ok<ASGLU>
k = n*4+1; %#ok<NASGU>

% pre-declare the variables the scripts assign that are read back afterwards,
% so the parser binds them as variables in this workspace
Ahat = []; Shat = []; KA = [];
ML = []; llike = []; lpri = []; lpost = []; store_lpost = []; den_psi = [];
store_A = []; store_Sig = []; store_nu = []; store_lam = []; store_pnu = [];
store_h = []; store_theta = []; store_psi = []; store_ppsi = [];
A_mean = []; Sig_mean = []; nu_mean = []; pnu_mean = []; h_mean = [];
theta_mean = []; CSV_std_mean = []; psi_mean = []; ppsi_mean = [];
countnu = []; counth = []; countrho = []; countpsi = [];
count_h = []; count_rho = []; count_psi = []; count_nu = [];

scripts = {'BVAR','BVAR_t','BVAR_CSV','BVAR_MA','BVAR_t_CSV','BVAR_t_MA', ...
    'BVAR_CSV_MA','BVAR_CSV_t_MA'};
resolved = which(scripts{model});
assert(strncmpi(resolved, tmp, numel(tmp)), ...
    '%s must resolve from the tempdir copy, got %s', scripts{model}, resolved);
resolved = which('construct_prior_A');
assert(strncmpi(resolved, tmp, numel(tmp)), ...
    'construct_prior_A must resolve from the tempdir copy, got %s', resolved);

rng(seed, 'twister');
switch model
    case 1, BVAR;
    case 2, BVAR_t;
    case 3, BVAR_CSV;
    case 4, BVAR_MA;
    case 5, BVAR_t_CSV;
    case 6, BVAR_t_MA;
    case 7, BVAR_CSV_MA;
    case 8, BVAR_CSV_t_MA;
end
s = rng;
close all   % the legacy scripts draw figure windows; drop them

names = {'Ahat','Shat','KA','ML','llike','lpri','lpost','store_lpost','den_psi', ...
    'store_A','store_Sig','store_nu','store_lam','store_pnu','store_h', ...
    'store_theta','store_psi','store_ppsi','A_mean','Sig_mean','nu_mean', ...
    'pnu_mean','h_mean','theta_mean','CSV_std_mean','psi_mean','ppsi_mean', ...
    'countnu','counth','countrho','countpsi','count_h','count_rho', ...
    'count_psi','count_nu'};
L = struct();
for kf = 1:numel(names)
    L.(names{kf}) = eval(names{kf});
end
L.rngstate = s.State;
end
