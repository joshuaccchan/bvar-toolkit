function test_oisv_equivalence
% seeded draw-for-draw equivalence of the functionized OISV full-sample
% pipeline (replications/chan_koop_yu2024_jbes_oisv/run_all.m +
% bvt.structural.b0_row_sampler / construct_Sigt + bvt.samplers.eq_svar_oi /
% eq_tri_cs / alp_tri_cs / horseshoe_kappa_psi + bvt.priors.* + bvt.sv.* +
% bvt.util.anormrnd) with the legacy workspace scripts SVARSV_MH ('OI', run at
% the default ordering) and CS_MH ('CS', run at the REVERSED ordering, so both
% rev_option branches are exercised), at small nsim from a tempdir copy.
% Asserts isequal on ALL stored draws (not just means), count_phi, the legacy
% script-tail posterior summaries, the six func_main_SVAR_v2 outputs
% (Mean_kappa/Median_kappa/Mean_Impact/Median_Impact/Sig_mean/Sig_median -
% the legacy side replicates func_main's post-processing with the tempdir copy
% of legacy construct_Sigt.m), and the terminal rng state (same rng call
% sequence).
%
% SOLE PATCH to the legacy scripts: SVARSV_MH.m line 24 carries the clock-seed
% line
%   randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
% ACTIVE (asserted: exactly one occurrence, not commented) - it re-seeds the
% global stream from the wall clock (irreproducible) AND switches MATLAB to
% the legacy v4/v5 generators, so it is removed from the tempdir copy.
% CS_MH.m carries the same line already COMMENTED OUT (line 49; asserted: its
% only occurrence is the commented one) and runs BYTE-VERBATIM. Every rng draw
% in both scripts sits after that point, so seeding rng(seed,'twister') just
% before dispatching aligns the whole run, chain-init draws included.
% func_main_SVAR_v2.m's setup (lines 5-31) is replicated in the harness below
% (run_legacy), which is also where the hard-coded nsim = 30000 / burnin =
% 5000 of main_SVAR_fullsample.m are overridden: the sampler scripts read
% nsim/burnin from the workspace, so the harness simply assigns the small
% values before dispatching the script.
root = getappdata(0, 'bvt_repo_root');
leg = fullfile(root, 'replications', 'chan_koop_yu2024_jbes_oisv', 'legacy');
repdir = fullfile(root, 'replications', 'chan_koop_yu2024_jbes_oisv');

nsim = 30; burnin = 10; seed = 20260902;    % ~0.9 s/sweep at T=706, n=20, p=13 - sized to keep the test under ~3 minutes

% --- tempdir with the patched OI script + byte-verbatim CS script/helper copies ---
tmp = tempname; mkdir(tmp);
ctmp = onCleanup(@() cleanup_tmp(tmp));   % rmpath BEFORE rmdir, warning-free
helpers = {'get_resid_var.m', 'get_C.m', 'getVbeta.m', 'vec.m', 'anormrnd.m', ...
    'sample_SV.m', 'sample_SV0para.m', 'sample_SVpara.m', 'construct_Sigt.m'};
for k = 1:numel(helpers)
    copyfile(fullfile(leg, 'utility', helpers{k}), fullfile(tmp, helpers{k}));
end
seedline = 'randn(''seed'',sum(clock*100)); rand(''seed'',sum(clock*1000));';

% SVARSV_MH.m: exactly one clock-seed occurrence, ACTIVE (not commented) - remove it
txt = fileread(fullfile(leg, 'SVARSV_MH.m'));
assert(numel(strfind(txt, seedline)) == 1, ...
    'expected exactly one clock-seed line in SVARSV_MH.m');
assert(isempty(strfind(txt, ['%' seedline])), ...
    'the SVARSV_MH.m clock-seed line is expected to be ACTIVE');
txt = strrep(txt, seedline, ...
    '% [clock-seed line removed by test_oisv_equivalence - the sole patch]');
fid = fopen(fullfile(tmp, 'SVARSV_MH.m'), 'w');
fwrite(fid, txt);
fclose(fid);

% CS_MH.m: no ACTIVE clock-seed line (its only occurrence is commented out) - byte-verbatim
txt = fileread(fullfile(leg, 'CS_MH.m'));
assert(numel(strfind(txt, ['%' seedline])) == 1 ...
    && numel(strfind(txt, seedline)) == 1, ...
    'expected CS_MH.m''s only clock-seed occurrence to be the commented-out one');
copyfile(fullfile(leg, 'CS_MH.m'), fullfile(tmp, 'CS_MH.m'));

addpath(tmp);   % removed by cleanup_tmp via ctmp, before the folder is deleted
addpath(repdir); cp2 = onCleanup(@() rmpath(repdir));

% both chan2021_ijf_mahp and this package define run_all - pin the resolution
resolved = which('run_all');
assert(strncmpi(resolved, repdir, numel(repdir)), ...
    'run_all must resolve from the OISV package, got %s', resolved);

% model, rev_option (flip), and the fields to compare
cases = { ...
    'OI', 0, {'store_kappa','store_hpara','store_alpha','store_B0','store_h', ...
              'count_phi','h_mean','B0_mean','B0_median','alpha_mean', ...
              'hpara_mean','kappa_mean','kappa_median','A_mean'}; ...
    'CS', 1, {'store_beta','store_alp','store_h','store_hpara','store_kappa', ...
              'count_phi','beta_mean','Beta_mean','alp_mean','alp_median', ...
              'h_mean','hpara_mean','kappa_mean','kappa_median'} };
funcmain_outs = {'Mean_kappa','Median_kappa','Mean_Impact','Median_Impact', ...
    'Sig_mean','Sig_median'};

for kc = 1:size(cases, 1)
    model = cases{kc, 1};
    rev = cases{kc, 2};
    if strcmp(model, 'OI')
        L = run_legacy_oi(rev, leg, tmp, nsim, burnin, seed);
    else
        L = run_legacy_cs(rev, leg, tmp, nsim, burnin, seed);
    end
    res = run_all(model, rev, nsim, burnin, seed);   % re-seeds itself from `seed`
    sC = rng;
    flds = cases{kc, 3};
    for kf = 1:numel(flds)
        assert(isequal(L.(flds{kf}), res.(flds{kf})), ...
            '%s: %s differs', model, flds{kf});
    end
    for kf = 1:numel(funcmain_outs)
        assert(isequal(L.(funcmain_outs{kf}), res.(funcmain_outs{kf})), ...
            '%s: %s differs', model, funcmain_outs{kf});
    end
    assert(isequal(L.rngstate, sC.State), '%s: rng call sequence differs', model);
end
end

function cleanup_tmp(tmp)
% deterministic teardown order: path entry first, then the folder itself
if any(strcmpi(strsplit(path, pathsep), tmp))
    rmpath(tmp);
end
if exist(tmp, 'dir')
    rmdir(tmp, 's');
end
end

function [Y0, Y, T, n, X, p] = prime_design(rev, leg)
% replicate func_main_SVAR_v2.m's setup (lines 8-31): p, data, variable
% ordering (rev_option), initial conditions, and the lagged design matrix
p = 13;
data = load(fullfile(leg, 'FRED_MD_20vars.csv'));       % legacy data, read-only
var_id1 = [4,6,12,13];
var_id2 = [1:3,5,7:11,14:20];
var_id = [var_id1, var_id2];
if rev == 1
    var_id = flip(var_id);
end
Y0 = data(1:24, var_id);
Y = data(25:end, var_id);
[T, n] = size(Y);
tmpY = [Y0(end-p+1:end, :); Y];
X = zeros(T, n*p);
for ii = 1:p
    X(:, (ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii, :);
end
X = [ones(T, 1) X];
end

function out = run_legacy_oi(rev, leg, tmp, nsim, burnin, seed) %#ok<INUSD> % nsim/burnin are read by the dispatched script from this workspace
% prime the workspace exactly as func_main_SVAR_v2.m does, dispatch the
% PATCHED SVARSV_MH from the tempdir, then replicate func_main's model-1
% post-processing (lines 36-46) with the tempdir LEGACY construct_Sigt copy
[Y0, Y, T, n, X, p] = prime_design(rev, leg);           %#ok<ASGLU> % Y0/Y/T/n/X/p are read by the dispatched script

% pre-declare the variables the script assigns that are read back afterwards,
% so the parser binds them as variables in this workspace
store_kappa = []; store_hpara = []; store_alpha = []; store_B0 = []; store_h = [];
count_phi = []; h_mean = []; B0_mean = []; B0_median = []; alpha_mean = [];
hpara_mean = []; kappa_mean = []; kappa_median = []; A_mean = [];

resolved = which('SVARSV_MH');
assert(strncmpi(resolved, tmp, numel(tmp)), ...
    'patched SVARSV_MH must resolve from the tempdir copy, got %s', resolved);
resolved = which('construct_Sigt');
assert(strncmpi(resolved, tmp, numel(tmp)), ...
    'construct_Sigt must resolve from the tempdir legacy copy, got %s', resolved);

rng(seed, 'twister');
SVARSV_MH;
s = rng;

    % func_main_SVAR_v2.m lines 36-46 (model 1), verbatim
Mean_kappa = kappa_mean; Median_kappa = median(store_kappa)';
Mean_Impact = reshape(B0_mean,n,n)';
B0_median = median(store_B0)'; Median_Impact = reshape(B0_median,n,n)';

store_Sig = zeros(T,n,n);
for isim = 1:nsim
    B0 = reshape(store_B0(isim,:),n,n)'; % stacked by rows
    h = squeeze(store_h(isim,:,:));
    store_Sig = store_Sig + construct_Sigt(h,B0);
end
Sig_mean = store_Sig/nsim; Sig_median = zeros(T,n,n);

out = struct('store_kappa', store_kappa, 'store_hpara', store_hpara, ...
    'store_alpha', store_alpha, 'store_B0', store_B0, 'store_h', store_h, ...
    'count_phi', count_phi, 'h_mean', h_mean, 'B0_mean', B0_mean, ...
    'B0_median', B0_median, 'alpha_mean', alpha_mean, 'hpara_mean', hpara_mean, ...
    'kappa_mean', kappa_mean, 'kappa_median', kappa_median, 'A_mean', A_mean, ...
    'Mean_kappa', Mean_kappa, 'Median_kappa', Median_kappa, ...
    'Mean_Impact', Mean_Impact, 'Median_Impact', Median_Impact, ...
    'Sig_mean', Sig_mean, 'Sig_median', Sig_median, 'rngstate', s.State);
end

function out = run_legacy_cs(rev, leg, tmp, nsim, burnin, seed) %#ok<INUSD> % nsim/burnin are read by the dispatched script from this workspace
% prime the workspace exactly as func_main_SVAR_v2.m does, dispatch the
% byte-verbatim CS_MH copy from the tempdir, then replicate func_main's
% model-2 post-processing (lines 49-62) with the tempdir LEGACY construct_Sigt
[Y0, Y, T, n, X, p] = prime_design(rev, leg);           %#ok<ASGLU> % Y0/Y/T/n/X/p are read by the dispatched script

store_beta = []; store_alp = []; store_h = []; store_hpara = []; store_kappa = [];
count_phi = []; beta_mean = []; Beta_mean = []; alp_mean = []; alp_median = [];
h_mean = []; hpara_mean = []; kappa_mean = []; kappa_median = [];

resolved = which('CS_MH');
assert(strncmpi(resolved, tmp, numel(tmp)), ...
    'CS_MH must resolve from the tempdir copy, got %s', resolved);

rng(seed, 'twister');
CS_MH;
s = rng;

    % func_main_SVAR_v2.m lines 49-62 (model 2), verbatim
Mean_kappa = kappa_mean; Median_kappa = median(store_kappa)';
A_id = nonzeros(tril(reshape(1:n^2,n,n),-1)'); A = eye(n); A(A_id) = alp_mean;
Mean_Impact = A;
A_id = nonzeros(tril(reshape(1:n^2,n,n),-1)'); A = eye(n); alp_median = median(store_alp)'; A(A_id) = alp_median;
Median_Impact = A;

store_Sig = zeros(T,n,n);
for isim = 1:nsim
    h = squeeze(store_h(isim,:,:));
    alp = squeeze(store_alp(isim,:))';
    A = eye(n); A(A_id) = alp;
    store_Sig = store_Sig + construct_Sigt(h,A);
end
Sig_mean = store_Sig/nsim; Sig_median = zeros(T,n,n);

out = struct('store_beta', store_beta, 'store_alp', store_alp, ...
    'store_h', store_h, 'store_hpara', store_hpara, 'store_kappa', store_kappa, ...
    'count_phi', count_phi, 'beta_mean', beta_mean, 'Beta_mean', Beta_mean, ...
    'alp_mean', alp_mean, 'alp_median', alp_median, 'h_mean', h_mean, ...
    'hpara_mean', hpara_mean, 'kappa_mean', kappa_mean, 'kappa_median', kappa_median, ...
    'Mean_kappa', Mean_kappa, 'Median_kappa', Median_kappa, ...
    'Mean_Impact', Mean_Impact, 'Median_Impact', Median_Impact, ...
    'Sig_mean', Sig_mean, 'Sig_median', Sig_median, 'rngstate', s.State);
end
