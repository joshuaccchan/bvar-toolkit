function test_mahp_equivalence
% seeded draw-for-draw equivalence of the functionized MAHP estimation pipeline
% (replications/chan2021_ijf_mahp/run_all.m + bvar.samplers.eq_gauss /
% gig_shrinkage / nu_psi_ng + bvar.priors.* + bvar.sv.ksc_rw_h0) with the legacy
% workspace scripts BVAR_MNG / BVAR_NG / BVAR_Minn, run at small nsim from a
% tempdir copy. Asserts isequal on ALL stored draws (not just means), on the
% posterior means, and on the terminal rng state (same rng call sequence).
%
% SOLE PATCH to the legacy scripts: the clock-seed line
%   randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
% (BVAR_MNG.m line 32, BVAR_NG.m line 31, BVAR_Minn.m line 25) is removed - it
% re-seeds the global stream from the wall clock (irreproducible) AND switches
% MATLAB to the legacy v4/v5 generators. Everything else runs byte-verbatim.
% main_BVAR.m's setup is replicated in the harness below (run_legacy), which is
% also where the hard-coded nsim = 10000 / burnin = 1000 are overridden: the
% sampler scripts read nsim/burnin from the workspace, so the harness simply
% assigns the small values before dispatching the script.
root = getappdata(0, 'bvar_repo_root');
leg = fullfile(root, 'replications', 'chan2021_ijf_mahp', 'legacy');
repdir = fullfile(root, 'replications', 'chan2021_ijf_mahp');

nsim = 200; burnin = 50; seed = 20260901;

% --- tempdir with the patched sampler scripts + byte-verbatim helper copies ---
tmp = tempname; mkdir(tmp);
ctmp = onCleanup(@() cleanup_tmp(tmp));   % rmpath BEFORE rmdir, warning-free
helpers = {'SVRW.m', 'getVtheta.m', 'get_C.m', 'get_resid_var.m', 'gigrnd.m', 'sample_nu_psi.m'};
for k = 1:numel(helpers)
    copyfile(fullfile(leg, helpers{k}), fullfile(tmp, helpers{k}));
end
seedline = 'randn(''seed'',sum(clock*100)); rand(''seed'',sum(clock*1000));';
scripts = {'BVAR_MNG.m', 'BVAR_NG.m', 'BVAR_Minn.m'};
for k = 1:numel(scripts)
    txt = fileread(fullfile(leg, scripts{k}));
    assert(numel(strfind(txt, seedline)) == 1, ...
        'expected exactly one clock-seed line in %s', scripts{k});
    txt = strrep(txt, seedline, ...
        '% [clock-seed line removed by test_mahp_equivalence - the sole patch]');
    fid = fopen(fullfile(tmp, scripts{k}), 'w');
    fwrite(fid, txt);
    fclose(fid);
end

addpath(tmp);   % removed by cleanup_tmp via ctmp, before the folder is deleted
addpath(repdir); cp2 = onCleanup(@() rmpath(repdir));

models = {'MNG', 'BVAR_MNG'; 'NG', 'BVAR_NG'; 'Minn', 'BVAR_Minn'};
for k = 1:size(models, 1)
    L = run_legacy(models{k, 2}, leg, tmp, nsim, burnin, seed);
    res = run_all(models{k, 1}, nsim, burnin, seed);   % re-seeds itself from `seed`
    sC = rng;
    assert(isequal(L.store_beta,  res.store_beta),  '%s: store_beta differs',  models{k, 1});
    assert(isequal(L.store_alp,   res.store_alp),   '%s: store_alp differs',   models{k, 1});
    assert(isequal(L.store_h,     res.store_h),     '%s: store_h differs',     models{k, 1});
    assert(isequal(L.store_Sigh,  res.store_Sigh),  '%s: store_Sigh differs',  models{k, 1});
    assert(isequal(L.store_kappa, res.store_kappa), '%s: store_kappa differs', models{k, 1});
    assert(isequal(L.beta_hat, res.beta_hat) && isequal(L.alp_hat, res.alp_hat) ...
        && isequal(L.h_hat, res.h_hat) && isequal(L.kappa_hat, res.kappa_hat), ...
        '%s: posterior means differ', models{k, 1});
    if ~strcmp(models{k, 1}, 'Minn')
        assert(isequal(L.store_nu_psi, res.store_nu_psi), '%s: store_nu_psi differs', models{k, 1});
        assert(isequal(L.count_nu_psi, res.count_nu_psi), '%s: count_nu_psi differs', models{k, 1});
    end
    assert(isequal(L.rngstate, sC.State), '%s: rng call sequence differs', models{k, 1});
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

function out = run_legacy(script_name, leg, tmp, nsim, burnin, seed) %#ok<INUSD> % nsim/burnin are read by the dispatched script from this workspace
% replicate main_BVAR.m's setup (lines 19-46) in this workspace - with the
% small nsim/burnin already assigned (the override point) - then dispatch the
% PATCHED sampler script exactly as main_BVAR.m lines 48-55 does.
p = 4;
raw = load(fullfile(leg, 'macrodata_Q_2018Q4.csv'));    % legacy data, read-only
data = raw(1:238, :);
var_id = [1,2,22,23,35,37,57,58,59,76,81,83,95,120,138,144,145,147,148,152,160,161,245];
Y0 = data(1:8, var_id);
Y = data(9:end, var_id);
[T, n] = size(Y);
tmpY = [Y0(end-p+1:end, :); Y];
Z = zeros(T, n*p);
for ii = 1:p
    Z(:, (ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii, :);
end
Z = [ones(T, 1) Z];                                     %#ok<NASGU> % read by the dispatched script
k_beta = n^2*p + n;                                     %#ok<NASGU>
k_alp = n*(n-1)/2;                                      %#ok<NASGU>
sig2 = get_resid_var(Y0, Y);                            % tempdir legacy copy
[C, idx_kappa1, idx_kappa2] = get_C(n, p, sig2);        %#ok<ASGLU> % tempdir legacy copy
ah = zeros(n, 1); Vh = 10*ones(n, 1);                   %#ok<NASGU>
nuh0 = 5*ones(n, 1); Sh0 = .01*ones(n, 1).*(nuh0-1);    %#ok<NASGU>
c01 = [1, 1/.04];                                       %#ok<NASGU>
c02 = [1, 1/.04^2];                                     %#ok<NASGU>
lam0_nu_psi = 1;                                        %#ok<NASGU>

% pre-declare the variables the script assigns that are read back afterwards,
% so the parser binds them as variables in this workspace
store_beta = []; store_alp = []; store_h = []; store_Sigh = []; store_kappa = [];
store_nu_psi = []; count_nu_psi = [];
beta_hat = []; alp_hat = []; h_hat = []; kappa_hat = [];

resolved = which(script_name);
assert(strncmpi(resolved, tmp, numel(tmp)), ...
    'patched %s must resolve from the tempdir copy, got %s', script_name, resolved);

rng(seed, 'twister');
switch script_name
    case 'BVAR_MNG',  BVAR_MNG;
    case 'BVAR_NG',   BVAR_NG;
    case 'BVAR_Minn', BVAR_Minn;
    otherwise, error('unknown legacy script %s', script_name);
end
s = rng;
close all   % the legacy script tails create histogram figures

out = struct('store_beta', store_beta, 'store_alp', store_alp, ...
    'store_h', store_h, 'store_Sigh', store_Sigh, 'store_kappa', store_kappa, ...
    'store_nu_psi', store_nu_psi, 'count_nu_psi', count_nu_psi, ...
    'beta_hat', beta_hat, 'alp_hat', alp_hat, 'h_hat', h_hat, ...
    'kappa_hat', kappa_hat, 'rngstate', s.State);
end
