function test_mlvarsv_equivalence
% seeded draw-for-draw equivalence of the functionized ml_varsv estimation
% pipeline (replications/chan2023_joe_mlvarsv/run_all.m + bvar.samplers.
% eq_var_redu_tri / alp_tri_cs / factor_fsv / eq_fsv_load + bvar.sv.svo_outlier /
% csv_armh / ksc_ar1_mean / sv_params / init_approx1N + bvar.priors.*) with the
% legacy workspace scripts VAR_NCP, VAR_CSV, VAR_ARSV_redu, VAR_FSV and
% VAR_ARSVO_redu, run from tempdir copies at small nsim on the legacy n = 7
% variable subset. Asserts isequal on ALL stored draws (not just means),
% counters, the script-tail summaries, and the terminal rng state.
%
% Sole patch to the legacy scripts: the four MCMC scripts carry
%   randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
% active (asserted: exactly one occurrence each, not commented) - it re-seeds
% from the wall clock and switches MATLAB to the v4/v5 generators, so it is
% removed. VAR_NCP.m has no such line (asserted) and runs byte-verbatim. Every
% rng draw sits after that point, so rng(seed,'twister') before dispatch aligns
% the whole run, chain init included.
%
% Scope: estimation only. cp_ml = 0 for models 2-5 so the utility/ml_var_*
% routines (a separate phase) are never entered; model 1 runs with cp_ml = 1
% because VAR_NCP.m's log-ML is inline and analytic.
root = getappdata(0, 'bvar_repo_root');
leg = fullfile(root, 'replications', 'chan2023_joe_mlvarsv', 'legacy');
repdir = fullfile(root, 'replications', 'chan2023_joe_mlvarsv');

nsim = 60; burnin = 20; seed = 20260903;
varid = [1,22,59,120,133,144,148];   % main_varsv.m line 34 (the commented n = 7 subset)

% --- tempdir: patched estimation scripts + byte-verbatim utility copies ---
tmp = tempname; mkdir(tmp);
ctmp = onCleanup(@() cleanup_tmp(tmp));   % rmpath BEFORE rmdir, warning-free
helpers = {'prior_Minn.m', 'prior_NCP.m', 'prior_B0.m', 'sample_SV.m', ...
    'sample_SVpara.m', 'sample_CSV.m', 'get_C.m', 'getARh_approx1N.m', ...
    'vec.m', 'vech.m', 'gigrnd.m', 'ldet.m', 'mgammaln.m'};
for k = 1:numel(helpers)
    copyfile(fullfile(leg, 'utility', helpers{k}), fullfile(tmp, helpers{k}));
end

seedline = 'randn(''seed'',sum(clock*100)); rand(''seed'',sum(clock*1000));';
patched = {'VAR_CSV.m', 'VAR_ARSV_redu.m', 'VAR_FSV.m', 'VAR_ARSVO_redu.m'};
for k = 1:numel(patched)
    txt = fileread(fullfile(leg, patched{k}));
    assert(numel(strfind(txt, seedline)) == 1, ...
        'expected exactly one clock-seed line in %s', patched{k});
    assert(isempty(strfind(txt, ['%' seedline])), ...
        'the %s clock-seed line is expected to be ACTIVE', patched{k});
    txt = strrep(txt, seedline, ...
        '% [clock-seed line removed by test_mlvarsv_equivalence - the sole patch]');
    fid = fopen(fullfile(tmp, patched{k}), 'w');
    fwrite(fid, txt);
    fclose(fid);
end
% VAR_NCP.m: no MCMC, no clock-seed line - byte-verbatim
txt = fileread(fullfile(leg, 'VAR_NCP.m'));
assert(isempty(strfind(txt, 'clock*100')), 'VAR_NCP.m is expected to carry no clock-seed line');
copyfile(fullfile(leg, 'VAR_NCP.m'), fullfile(tmp, 'VAR_NCP.m'));

addpath(tmp);   % removed by cleanup_tmp via ctmp, before the folder is deleted
addpath(repdir); cp2 = onCleanup(@() rmpath(repdir)); %#ok<NASGU>

% four packages define run_all - pin the resolution
resolved = which('run_all');
assert(strncmpi(resolved, repdir, numel(repdir)), ...
    'run_all must resolve from the ml_varsv package, got %s', resolved);

% run_all's unqualified gigrnd resolves to the tempdir LEGACY copy here, not to
% third_party/gigrnd.m; assert the two are code-identical so the guarantee does
% not depend silently on that (see the step-5 path note in tests/variant_map.md)
assert(isequal(code_of(fullfile(tmp, 'gigrnd.m')), ...
               code_of(fullfile(root, 'third_party', 'gigrnd.m'))), ...
    'legacy gigrnd.m and third_party/gigrnd.m differ in code, not just comments');

% VAR_CSV.m ends with a figure - keep it off screen and clean up
vis = get(groot, 'defaultFigureVisible');
set(groot, 'defaultFigureVisible', 'off');
cfig = onCleanup(@() restore_figs(vis)); %#ok<NASGU>

% fields compared per model
f{1} = {'A_tilde','K_A','A_hat','S_hat','lml'};
f{2} = {'store_Sig','store_a','store_h','store_hpara','store_kappa', ...
        'count_h','count_phi','h_mean','hpara_mean','kappa_mean','CSV_std_mean'};
f{3} = {'store_alp','store_beta','store_h','store_hpara','store_kappa', ...
        'count_phi','alp_mean','A_mean','beta_mean','h_mean','hpara_mean','kappa_mean'};
f{4} = {'store_l','store_A','store_h','store_F','store_hpara','store_kappa', ...
        'count_phi','l_mean','L_mean','A_mean','h_mean','F_mean','hpara_mean','kappa_mean'};
f{5} = {'store_alp','store_beta','store_h','store_o','store_po', ...
        'store_hpara','store_kappa','count_phi','alp_mean','A_mean', ...
        'beta_mean','h_mean','hpara_mean','kappa_mean','o_mean','po_mean'};

% model, is_kappafixed, is_kappasym, varid. Every model runs at both switch
% settings it reads (VAR-NCP reads neither; VAR-CSV reads only is_kappafixed -
% main_varsv.m 68-73 has no is_kappasym branch for it). The last three rows
% repeat the default configuration at the paper's active n = 15 selection, so
% the per-equation index arithmetic is checked at a second dimension.
varid15 = [1,22,59,120,133,144,148,2,35,57,81,95,152,160,245];   % main_varsv.m line 35
cases = { ...
  1, false, false, varid; ...
  2, false, false, varid; ...
  2, true,  false, varid; ...
  3, false, false, varid; ...
  3, true,  false, varid; ...
  3, false, true,  varid; ...
  4, false, false, varid; ...
  4, true,  false, varid; ...
  4, false, true,  varid; ...
  5, false, false, varid; ...
  5, true,  false, varid; ...
  5, false, true,  varid; ...
  3, false, false, varid15; ...
  4, false, false, varid15; ...
  5, false, false, varid15 };

saw_outlier = false;
for kc = 1:size(cases, 1)
    imodel = cases{kc,1}; kfix = cases{kc,2}; ksym = cases{kc,3}; vid = cases{kc,4};
    tag = sprintf('model %d (n=%d, kappafixed=%d, kappasym=%d)', ...
        imodel, numel(vid), kfix, ksym);

    L = run_legacy(imodel, kfix, ksym, leg, tmp, vid, nsim, burnin, seed);
    res = run_all(imodel, kfix, ksym, nsim, burnin, seed, vid);   % re-seeds itself
    sC = rng;

    flds = f{imodel};
    for kf = 1:numel(flds)
        assert(isequal(L.(flds{kf}), res.(flds{kf})), '%s: %s differs', tag, flds{kf});
    end
    assert(isequal(L.rngstate, sC.State), '%s: rng call sequence differs', tag);
    assert(strcmp(L.model_name, res.model_name), '%s: model_name differs', tag);

    if imodel == 5
        saw_outlier = saw_outlier || any(res.store_o(:) > 1);
    end
end

% the SVO runs must actually visit outlier states, or the o block is untested
assert(saw_outlier, 'SVO: no outlier ever drawn - the o block is not exercised');
end

% -------------------------------------------------------------------------
function restore_figs(vis)
close all force
set(groot, 'defaultFigureVisible', vis);
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

function c = code_of(f)
% file contents with whole-line comments and blank lines removed
lines = strsplit(fileread(f), newline);
keep = ~cellfun(@(s) isempty(strtrim(s)) || strncmp(strtrim(s), '%', 1), lines);
c = strjoin(lines(keep), newline);
end

% -------------------------------------------------------------------------
function out = run_legacy(imodel, is_kappafixed, is_kappasym, leg, tmp, varid, nsim, burnin, seed) %#ok<INUSD>
% replicate main_varsv.m lines 33-137 in this workspace and dispatch the
% (patched) legacy script from the tempdir. nsim/burnin/is_kappa* are read by
% the dispatched script from here, hence the INUSD suppression above.
p = 4;          % main_varsv.m 24
r = 2;          % 25
    % VAR_NCP's ML is inline and analytic; models 2-5 would enter
    % utility/ml_var_* (a separate phase), so their cp_ml stays off
cp_ml = (imodel == 1);

    % load data [33-53]
data_all = load(fullfile(leg, 'macrodata_Q_2019Q4.csv'));
data = data_all(:,varid);
Y0 = data(1:8,:);
Y = data(9:end,:);
[T,n] = size(Y);            %#ok<ASGLU>
y = reshape(Y',n*T,1);      %#ok<NASGU>
tmpY = [Y0(end-p+1:end,:); Y];
X = zeros(T,n*p);
for i=1:p
    X(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end
X = [ones(T,1) X];
k_beta = n*(n-1)/2;         %#ok<NASGU>
k_alp = n^2*p + n;
k = k_alp/n;                %#ok<NASGU>

    % priors [56-122]
kappa = .2^2;
kappa1 = .2^2; kappa3 = 100; kappa4 = .2^2;
if is_kappasym
    kappa2 = kappa1;        %#ok<NASGU>
else
    kappa2 = (.2^2)^2;      %#ok<NASGU>
end
sig2_hat = []; U_hat = []; Hyper = struct();  %#ok<NASGU>
switch imodel
    case 1
        model_name = 'VAR-NCP';
        [Hyper.A0,Hyper.VA,Hyper.nu0,Hyper.S0] = prior_NCP(p,kappa,kappa3,Y0,Y);
    case 2
        if is_kappafixed
            model_name = 'VAR-CSV-fixed-kappa';
        else
            model_name = 'VAR-CSV';
        end
        [Hyper.A0,Hyper.VA,Hyper.nu0,Hyper.S0] = prior_NCP(p,kappa,kappa3,Y0,Y);
        Hyper.c0 = [1,1/.2^2];
        Hyper.nuh = 3; Hyper.Sh = .1*(Hyper.nuh-1);
        Hyper.phi0 = .98; Hyper.Vphi = .05^2;
    case 3
        if is_kappafixed
            model_name = 'VAR-SV-fixed-kappa';
        elseif is_kappasym
            model_name = 'VAR-SV-sym';
        else
            model_name = 'VAR-SV';
        end
        [Hyper.alp0,Hyper.Valp,sig2_hat,U_hat] = prior_Minn(p,kappa1,kappa2,kappa3,Y0,Y); %#ok<ASGLU>
        [Hyper.beta0,Hyper.Vbeta] = prior_B0(Y0,Y,kappa4);
        Hyper.c0 = [1,1/.2^2; 1,1/.2^2; 1,1];
        Hyper.nuh = 3*ones(n,1); Hyper.Sh = .1*(Hyper.nuh-1);
        Hyper.mu0 = zeros(n,1); Hyper.Vmu = 100*ones(n,1);
        Hyper.phi0 = .98*ones(n,1); Hyper.Vphi = .05^2*ones(n,1);
    case 4
        if is_kappafixed
            model_name = 'VAR-FSV-fixed-kappa';
        elseif is_kappasym
            model_name = 'VAR-FSV-sym';
        else
            model_name = 'VAR-FSV';
        end
        [Hyper.alp0,Hyper.Valp,sig2_hat] = prior_Minn(p,kappa1,kappa2,kappa3,Y0,Y);
        Hyper.c0 = [1,1/.2^2; 1,1/.2^2;];
        Hyper.nuh = 3*ones(n+r,1); Hyper.Sh = .1*(Hyper.nuh-1);
        Hyper.mu0 = zeros(n+r,1); Hyper.Vmu = 100*ones(n+r,1);
        Hyper.phi0 = .98*ones(n+r,1); Hyper.Vphi = .05^2*ones(n+r,1);
        Hyper.l0 = 0;
        Hyper.Vl = 1;
    case 5
        if is_kappafixed
            model_name = 'VAR-SVO-fixed-kappa';
        elseif is_kappasym
            model_name = 'VAR-SVO-sym';
        else
            model_name = 'VAR-SVO';
        end
        [Hyper.alp0,Hyper.Valp,sig2_hat,U_hat] = prior_Minn(p,kappa1,kappa2,kappa3,Y0,Y); %#ok<ASGLU>
        [Hyper.beta0,Hyper.Vbeta] = prior_B0(Y0,Y,kappa4);
        Hyper.c0 = [1,1/.2^2; 1,1/.2^2; 1,1];
        Hyper.nuh = 3*ones(n,1); Hyper.Sh = .1*(Hyper.nuh-1);
        Hyper.mu0 = zeros(n,1); Hyper.Vmu = 100*ones(n,1);
        Hyper.phi0 = .98*ones(n,1); Hyper.Vphi = .05^2*ones(n,1);
        Hyper.p0a = 10/4; Hyper.p0b = (1-1/16)*40;
end

% pre-declare the variables the dispatched script assigns and we read back,
% so the parser binds them as variables in this workspace
A_tilde = []; K_A = []; A_hat = []; S_hat = []; lml = [];
store_Sig = []; store_a = []; store_h = []; store_hpara = []; store_kappa = [];
store_alp = []; store_beta = []; store_l = []; store_A = []; store_F = [];
store_o = []; store_po = [];
count_h = []; count_phi = [];
h_mean = []; hpara_mean = []; kappa_mean = []; CSV_std_mean = [];
alp_mean = []; A_mean = []; beta_mean = []; l_mean = []; L_mean = []; F_mean = [];
o_mean = []; po_mean = [];

scripts = {'VAR_NCP', 'VAR_CSV', 'VAR_ARSV_redu', 'VAR_FSV', 'VAR_ARSVO_redu'};
resolved = which(scripts{imodel});
assert(strncmpi(resolved, tmp, numel(tmp)), ...
    '%s must resolve from the tempdir copy, got %s', scripts{imodel}, resolved);

rng(seed, 'twister');
switch imodel
    case 1, VAR_NCP;
    case 2, VAR_CSV;
    case 3, VAR_ARSV_redu;
    case 4, VAR_FSV;
    case 5, VAR_ARSVO_redu;
end
s = rng;

out = struct('A_tilde',A_tilde, 'K_A',K_A, 'A_hat',A_hat, 'S_hat',S_hat, 'lml',lml, ...
    'store_Sig',store_Sig, 'store_a',store_a, 'store_h',store_h, ...
    'store_hpara',store_hpara, 'store_kappa',store_kappa, 'store_alp',store_alp, ...
    'store_beta',store_beta, 'store_l',store_l, 'store_A',store_A, ...
    'store_F',store_F, 'store_o',store_o, 'store_po',store_po, ...
    'count_h',count_h, 'count_phi',count_phi, 'h_mean',h_mean, ...
    'hpara_mean',hpara_mean, 'kappa_mean',kappa_mean, 'CSV_std_mean',CSV_std_mean, ...
    'alp_mean',alp_mean, 'A_mean',A_mean, 'beta_mean',beta_mean, ...
    'l_mean',l_mean, 'L_mean',L_mean, 'F_mean',F_mean, ...
    'o_mean',o_mean, 'po_mean',po_mean, ...
    'model_name',model_name, 'rngstate',s.State);
end
