function test_hybtvp_equivalence
% seeded draw-for-draw equivalence of the functionized hybrid TVP-VAR pipeline
% (replications/chan2023_jbes_hybtvp/run_all.m + bvar.samplers.eq_hyb_tvp +
% bvar.sv.ksc_rw_h0 + bvar.priors.vtheta / minnesota_C / resid_var_allvars_ridge
% + bvar.util.gam_mode) with the legacy script main_HYB_TVPSV.m, run from a
% tempdir copy at small nsim. Asserts isequal on ALL stored draws, on the
% post-processing the legacy script prints (log Bayes factors, gamma summaries),
% and on the terminal rng state.
%
% PATCHES TO THE LEGACY SCRIPT. Unlike the other packages this one is a single
% self-contained script that opens with `clear`, so its run length cannot be
% injected through the caller's workspace. Four edits are made to the tempdir
% copy, each asserted to match exactly once:
%   1. the clock-seed line (main 85), removed - it re-seeds from the wall clock
%      and switches MATLAB to the v4/v5 generators;
%   2. nsim (main 16), reduced from 50000;
%   3. burnin (main 17), reduced from 1000;
%   4. `clear; clc;` (main 14), removed - run from a function it would wipe the
%      harness's own working-directory bookkeeping mid-run. It computes nothing.
% None of the four touches an arithmetic line.
% The var_id line is NOT patched. Instead the n = 3 selection on the commented
% line 22 is used by patching which line is active, so the test exercises the
% real index arithmetic at a dimension small enough to run in seconds.
%
% Note the seeding order this checks. main_HYB_TVPSV.m runs its `initialize`
% block - which itself draws n volatility paths - BEFORE the clock-seed line.
% With that line removed, one continuous stream covers initialization and the
% loop, and run_all reproduces exactly that: rng(seed) once, then initialize,
% then the sweeps.
root = getappdata(0, 'bvar_repo_root');
leg = fullfile(root, 'replications', 'chan2023_jbes_hybtvp', 'legacy');
repdir = fullfile(root, 'replications', 'chan2023_jbes_hybtvp');

nsim = 25; burnin = 10; seed = 20260905;

% --- tempdir: patched main script, verbatim utilities, and the data file ---
tmp = tempname; mkdir(tmp); mkdir(fullfile(tmp, 'utility'));
ctmp = onCleanup(@() cleanup_tmp(tmp));
u = dir(fullfile(leg, 'utility', '*.m'));
for k = 1:numel(u)
    copyfile(fullfile(leg, 'utility', u(k).name), fullfile(tmp, 'utility', u(k).name));
end
copyfile(fullfile(leg, 'macrodata_Q_2018Q4.csv'), fullfile(tmp, 'macrodata_Q_2018Q4.csv'));

txt = fileread(fullfile(leg, 'main_HYB_TVPSV.m'));
seedline = 'randn(''seed'',sum(clock*100)); rand(''seed'',sum(clock*1000));';
txt = patch_once(txt, seedline, ...
    '% [clock-seed line removed by test_hybtvp_equivalence]', 'clock-seed line');
txt = patch_once(txt, 'clear; clc;', ...
    '% [clear removed by test_hybtvp_equivalence]', 'clear statement');
txt = patch_once(txt, 'nsim = 50000; ', sprintf('nsim = %d;', nsim), 'nsim');
txt = patch_once(txt, 'burnin = 1000;', sprintf('burnin = %d;', burnin), 'burnin');
% activate the n = 3 selection (line 22) and retire the n = 6 one (line 23)
txt = patch_once(txt, '% var_id = [1,95,59]; % n = 3', 'var_id = [1,95,59];', 'n=3 selection');
txt = patch_once(txt, 'var_id = [1,95,59,144,22,133];   % n = 6', ...
    '% n = 6 selection retired by the test', 'n=6 selection');
fid = fopen(fullfile(tmp, 'main_HYB_TVPSV.m'), 'w'); fwrite(fid, txt); fclose(fid);

addpath(repdir); cp2 = onCleanup(@() rmpath(repdir)); %#ok<NASGU>

% four packages define run_all - pin the resolution
resolved = which('run_all');
assert(strncmpi(resolved, repdir, numel(repdir)), ...
    'run_all must resolve from the hybtvp package, got %s', resolved);

L = run_legacy(tmp, seed);
res = run_all(nsim, burnin, seed, [1,95,59], false, false);
sC = rng;

flds = {'store_alp','store_beta','store_h','store_Sigbeta','store_Sigalp', ...
        'store_Sigh','store_beta0','store_alp0','store_h0','store_p0', ...
        'store_gam','store_kappa','store_lpostgam'};
for k = 1:numel(flds)
    assert(isequal(L.(flds{k}), res.(flds{k})), '%s differs', flds{k});
end
assert(isequal(L.lBF, res.lBF), 'log Bayes factors differ');
assert(isequal(L.gam_hat, res.gam_hat), 'gam_hat differs');
assert(isequal(L.gam_mode, res.gam_mode), 'gam_mode differs');
assert(isequal(L.rngstate, sC.State), 'rng call sequence differs');

% the run must actually exercise time variation, or eq_hyb_tvp's interesting
% branches are never entered and the comparison proves little
assert(any(res.store_gam(:) == 1), 'no equation ever drew a time-varying block');
assert(any(res.store_gam(:) == 0), 'every draw was time-varying - the (0,0) branch is untested');
end

% -------------------------------------------------------------------------
function txt = patch_once(txt, from, to, what)
n = numel(strfind(txt, from));
assert(n == 1, 'expected exactly one occurrence of the %s, found %d', what, n);
txt = strrep(txt, from, to);
end

function cleanup_tmp(tmp)
% The legacy script's addpath('./utility') leaves an entry behind; drop every
% path entry under tmp before removing the folder.
entries = strsplit(path, pathsep);
for k = 1:numel(entries)
    if strncmpi(entries{k}, tmp, numel(tmp))
        rmpath(entries{k});
    end
end
if exist(tmp, 'dir')
    rmdir(tmp, 's');
end
end

% -------------------------------------------------------------------------
function out = run_legacy(tmpdir_, seed_)
% Run the patched legacy script from the tempdir as the working directory - its
% addpath('./utility') and its bare data filename are both relative. Locals here
% carry a trailing underscore so they cannot collide with the script's own
% variables (it assigns `tmp`, among others, in its printing loop).
od_ = cd(tmpdir_);

% Pre-declare everything read back below: the parser must bind these as
% variables, since the script assigns them through evalc.
store_alp = []; store_beta = []; store_h = []; store_Sigbeta = [];
store_Sigalp = []; store_Sigh = []; store_beta0 = []; store_alp0 = [];
store_h0 = []; store_p0 = []; store_gam = []; store_kappa = [];
store_lpostgam = []; lBF = []; gam_hat = []; gam_mode = [];

resolved_ = which('main_HYB_TVPSV');
assert(strncmpi(resolved_, tmpdir_, numel(tmpdir_)), ...
    'main_HYB_TVPSV must resolve from the tempdir copy, got %s', resolved_);

rng(seed_, 'twister');
try
    evalc('main_HYB_TVPSV');   % suppress its progress and summary printing
catch err_
    cd(od_);
    rethrow(err_);
end
s_ = rng;
cd(od_);

out = struct('store_alp',store_alp, 'store_beta',store_beta, 'store_h',store_h, ...
    'store_Sigbeta',store_Sigbeta, 'store_Sigalp',store_Sigalp, ...
    'store_Sigh',store_Sigh, 'store_beta0',store_beta0, 'store_alp0',store_alp0, ...
    'store_h0',store_h0, 'store_p0',store_p0, 'store_gam',store_gam, ...
    'store_kappa',store_kappa, 'store_lpostgam',store_lpostgam, ...
    'lBF',lBF, 'gam_hat',gam_hat, 'gam_mode',gam_mode, 'rngstate',s_.State);
end
