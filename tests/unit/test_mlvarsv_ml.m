function test_mlvarsv_ml
% seeded bitwise equivalence of the four extracted ml_varsv marginal-likelihood
% routines (bvar.ml.mlvarsv_csv / mlvarsv_arsv_redu / mlvarsv_fsv /
% mlvarsv_arsvo_redu, on top of bvar.ml.isden_arss and the three density
% utilities) with the legacy utility/ml_var_*.m routines, run from tempdir
% copies at small nsim and M on the full-length sample.
%
% Each configuration runs the legacy estimation script with cp_ml = 1 - so the
% legacy ML routine executes inside it, on the estimation's own rng stream -
% and the functionized pipeline (run_all + the core ML function) from the same
% seed, then asserts isequal on the stored draws, on lml and lmlstd, on
% store_w where the legacy script keeps it, and on the terminal rng state.
%
% VAR-SVO is compared under 'bugcompat', true: its legacy ML carries three
% defects in the outlier block (prior mass over 32 atoms instead of 31, the
% linear indexing of o_hat, and the missing -n*sum(log(o)) Jacobian). The
% corrected default is additionally asserted to consume the identical rng
% stream and to move the weights by exactly the three corrections.
%
% Sole patch to the legacy scripts: the four MCMC scripts' clock-seed line, as
% in test_mlvarsv_equivalence. The ml_var_* routines and every utility run
% byte-verbatim. Model 1 has no ML routine (its log-ML is inline and analytic
% in VAR_NCP.m) and is covered by test_mlvarsv_equivalence.
%
% T = 234 (the full sample) is deliberate: the o_hat linear-index defect only
% has its published shape when T >= 32.
root = getappdata(0, 'bvar_repo_root');
leg = fullfile(root, 'replications', 'chan2023_joe_mlvarsv', 'legacy');
repdir = fullfile(root, 'replications', 'chan2023_joe_mlvarsv');

nsim = 40; burnin = 10; seed = 20260903;
M = 100;            % rounded to a multiple of 50 inside every routine
flag_marg = 2;      % main_varsv.m line 30
varid = [1,22,59,120];   % n = 4: the first four of the commented n = 7 selection

% --- tempdir: patched estimation scripts + byte-verbatim utility copies ---
tmp = tempname; mkdir(tmp);
ctmp = onCleanup(@() cleanup_tmp(tmp));
helpers = {'prior_Minn.m', 'prior_NCP.m', 'prior_B0.m', 'sample_SV.m', ...
    'sample_SVpara.m', 'sample_CSV.m', 'get_C.m', 'getARh_approx1N.m', ...
    'vec.m', 'vech.m', 'gigrnd.m', 'ldet.m', 'mgammaln.m', 'SURform2.m', ...
    'tnormrnd.m', 'getISden_ARSS.m', 'lgampdf.m', 'ltnormpdf.m', 'lmvnpdf_pcn.m', ...
    'ml_var_csv.m', 'ml_var_arsv_redu.m', 'ml_var_fsv.m', 'ml_var_arsvo_redu.m'};
for k = 1:numel(helpers)
    copyfile(fullfile(leg, 'utility', helpers{k}), fullfile(tmp, helpers{k}));
end

seedline = 'randn(''seed'',sum(clock*100)); rand(''seed'',sum(clock*1000));';
patched = {'VAR_CSV.m', 'VAR_ARSV_redu.m', 'VAR_FSV.m', 'VAR_ARSVO_redu.m'};
for k = 1:numel(patched)
    txt = fileread(fullfile(leg, patched{k}));
    assert(numel(strfind(txt, seedline)) == 1, ...
        'expected exactly one clock-seed line in %s', patched{k});
    txt = strrep(txt, seedline, '% [clock-seed line removed by test_mlvarsv_ml - the sole patch]');
    fid = fopen(fullfile(tmp, patched{k}), 'w');
    fwrite(fid, txt);
    fclose(fid);
end
% the ML routines must be seed-line free (they run byte-verbatim)
mlfiles = {'ml_var_csv.m', 'ml_var_arsv_redu.m', 'ml_var_fsv.m', 'ml_var_arsvo_redu.m'};
for k = 1:numel(mlfiles)
    txt = fileread(fullfile(leg, 'utility', mlfiles{k}));
    assert(~contains(txt, 'clock*100'), '%s is expected to carry no clock-seed line', mlfiles{k});
    assert(any(strcmp(mlfiles{k}, helpers)), '%s must be copied to the tempdir', mlfiles{k});
end

addpath(tmp);
addpath(repdir); cp2 = onCleanup(@() rmpath(repdir)); %#ok<NASGU>

resolved = which('run_all');
assert(strncmpi(resolved, repdir, numel(repdir)), ...
    'run_all must resolve from the ml_varsv package, got %s', resolved);
assert(isequal(code_of(fullfile(tmp, 'gigrnd.m')), ...
               code_of(fullfile(root, 'third_party', 'gigrnd.m'))), ...
    'legacy gigrnd.m and third_party/gigrnd.m differ in code, not just comments');

% VAR_CSV.m ends with a figure
vis = get(groot, 'defaultFigureVisible');
set(groot, 'defaultFigureVisible', 'off');
cfig = onCleanup(@() restore_figs(vis)); %#ok<NASGU>

stores{2} = {'store_Sig','store_a','store_h','store_hpara','store_kappa'};
stores{3} = {'store_alp','store_beta','store_h','store_hpara','store_kappa'};
stores{4} = {'store_l','store_A','store_h','store_F','store_hpara','store_kappa'};
stores{5} = {'store_alp','store_beta','store_h','store_o','store_po','store_hpara','store_kappa'};

% model, is_kappafixed, is_kappasym. Every branch of every routine's
% cprior/gIS switch is visited (VAR-CSV reads only is_kappafixed).
cases = { 2, false, false; ...
          2, true,  false; ...
          3, false, false; ...
          3, true,  false; ...
          3, false, true;  ...
          4, false, false; ...
          4, true,  false; ...
          4, false, true;  ...
          5, false, false; ...
          5, true,  false; ...
          5, false, true };

for kc = 1:size(cases,1)
    imodel = cases{kc,1}; kfix = cases{kc,2}; ksym = cases{kc,3};
    tag = sprintf('model %d (kappafixed=%d, kappasym=%d)', imodel, kfix, ksym);

    L = run_legacy(imodel, kfix, ksym, leg, tmp, varid, nsim, burnin, seed, M, flag_marg);

    est = run_all(imodel, kfix, ksym, nsim, burnin, seed, varid);   % re-seeds itself
    s_est = rng;
    [lml, lmlstd, ml] = core_ml(imodel, est, M, flag_marg, true);
    s_end = rng;

    flds = stores{imodel};
    for kf = 1:numel(flds)
        assert(isequal(L.(flds{kf}), est.(flds{kf})), '%s: %s differs', tag, flds{kf});
    end
    assert(isequal(L.lml, lml), '%s: lml differs (%.15g vs %.15g)', tag, L.lml, lml);
    assert(isequal(L.lmlstd, lmlstd), '%s: lmlstd differs', tag);
    assert(isequal(L.rngstate, s_end.State), '%s: rng call sequence differs', tag);
    if imodel == 4    % the only script that keeps the third output
        assert(isequal(L.store_w, ml.store_w), '%s: store_w differs', tag);
    end
    assert(isfinite(lml) && isfinite(lmlstd) && lmlstd > 0, '%s: degenerate ML output', tag);

    if imodel == 3 && ~kfix && ~ksym
        % the prior blocks are rebuilt inside from each kappa draw, so the
        % alp0/Valp/beta0/Vbeta the caller passes cannot matter (the legacy
        % hands over the final sweep's, run_all the untouched originals)
        H2 = est.Hyper;
        H2.Valp = 3*H2.Valp; H2.alp0 = H2.alp0 + 1;
        H2.Vbeta = 7*H2.Vbeta; H2.beta0 = H2.beta0 + 1;
        est2 = est; est2.Hyper = H2;
        rng(s_est);
        lml2 = core_ml(imodel, est2, M, flag_marg, true);
        assert(isequal(lml2, lml), '%s: the passed-in prior blocks must not matter', tag);
    end

    if imodel == 4 && ~kfix && ~ksym
        % flag_marg = 1 (sig2 kept as a drawn parameter) is unreachable from
        % main_varsv.m but implemented in ml_var_fsv.m alone - check it directly
        rng(s_est);
        [l1,s1,w1] = ml_var_fsv(est.X,est.Y,est.Y0,M,est.Hyper,1, ...
            est.store_h,est.store_hpara,est.store_l,est.store_kappa,kfix,ksym);
        sa = rng;
        rng(s_est);
        [l2,s2,d2] = bvar.ml.mlvarsv_fsv(est.X,est.Y,est.Y0,M,est.Hyper,1, ...
            est.store_h,est.store_hpara,est.store_l,est.store_kappa,kfix,ksym);
        sb = rng;
        assert(isequal(l1,l2) && isequal(s1,s2) && isequal(w1,d2.store_w) ...
            && isequal(sa.State,sb.State), '%s: flag_marg = 1 differs', tag);
        assert(~isequal(l1,lml), '%s: flag_marg 1 and 2 must differ', tag);
    end

    if imodel == 5
        % corrected default: same stream, weights moved by exactly the three
        % corrections (the c1 term is re-associated, hence the tolerance)
        rng(s_est);
        [lmlc, ~, mlc] = core_ml(imodel, est, M, flag_marg, false);
        s_cor = rng;
        assert(isequal(s_cor.State, s_end.State), '%s: the corrections must consume no rng', tag);
        assert(ml.n_prior_atoms == 32 && mlc.n_prior_atoms == 31, '%s: prior atom counts', tag);
        assert(all(ml.store_lJ_o == 0), '%s: bugcompat must apply no o Jacobian', tag);
        assert(any(mlc.store_lJ_o < 0), '%s: no outlier drawn - the o block is not exercised', tag);
        pred = mlc.store_lJ_o + (mlc.store_lr_o - ml.store_lr_o);
        assert(max(abs((mlc.store_w - ml.store_w) - pred)) < 1e-6, ...
            '%s: the corrected-vs-legacy weight gap is not the three corrections', tag);
        assert(~isequal(lmlc, lml), '%s: the corrected ML must differ', tag);
    end
end

% run_ml wires estimation and ML together on one stream: same numbers again
res = run_ml(2, false, false, nsim, burnin, seed, varid, 'M', M, 'flag_marg', flag_marg);
L2 = run_legacy(2, false, false, leg, tmp, varid, nsim, burnin, seed, M, flag_marg);
assert(isequal(res.lml, L2.lml) && isequal(res.lmlstd, L2.lmlstd), 'run_ml: lml differs');
end

% -------------------------------------------------------------------------
function [lml, lmlstd, ml] = core_ml(imodel, est, M, flag_marg, bugcompat)
Y = est.Y; X = est.X; Y0 = est.Y0; H = est.Hyper;
kfix = est.is_kappafixed; ksym = est.is_kappasym;
switch imodel
    case 2
        [lml,lmlstd,ml] = bvar.ml.mlvarsv_csv(X,Y,Y0,M,H, ...
            est.store_h,est.store_hpara,est.store_kappa,kfix);
    case 3
        [lml,lmlstd,ml] = bvar.ml.mlvarsv_arsv_redu(X,Y,Y0,M,H,flag_marg, ...
            est.store_h,est.store_beta,est.store_hpara,est.store_kappa,kfix,ksym);
    case 4
        [lml,lmlstd,ml] = bvar.ml.mlvarsv_fsv(X,Y,Y0,M,H,flag_marg, ...
            est.store_h,est.store_hpara,est.store_l,est.store_kappa,kfix,ksym);
    case 5
        [lml,lmlstd,ml] = bvar.ml.mlvarsv_arsvo_redu(X,Y,Y0,M,H,flag_marg, ...
            est.store_h,est.store_beta,est.store_hpara,est.store_kappa, ...
            est.store_o,est.store_po,est.o_grid,kfix,ksym,'bugcompat',bugcompat);
end
end

% -------------------------------------------------------------------------
function restore_figs(vis)
close all force
set(groot, 'defaultFigureVisible', vis);
end

function cleanup_tmp(tmp)
if any(strcmpi(strsplit(path, pathsep), tmp))
    rmpath(tmp);
end
if exist(tmp, 'dir')
    rmdir(tmp, 's');
end
end

function c = code_of(f)
lines = strsplit(fileread(f), newline);
keep = ~cellfun(@(s) isempty(strtrim(s)) || strncmp(strtrim(s), '%', 1), lines);
c = strjoin(lines(keep), newline);
end

% -------------------------------------------------------------------------
function out = run_legacy(imodel, is_kappafixed, is_kappasym, leg, tmp, varid, nsim, burnin, seed, M, flag_marg) %#ok<INUSD>
% main_varsv.m lines 33-137 in this workspace, with cp_ml = 1, dispatching the
% (patched) legacy script from the tempdir; nsim/burnin/M/flag_marg/cp_ml and
% the is_kappa* switches are read by the dispatched script from here.
p = 4;          % main_varsv.m 24
r = 2;          % 25
cp_ml = 1;      % 29 - the ML routine runs inside the script

data_all = load(fullfile(leg, 'macrodata_Q_2019Q4.csv'));
data = data_all(:,varid);
Y0 = data(1:8,:);
Y = data(9:end,:);
[T,n] = size(Y);            %#ok<ASGLU>
tmpY = [Y0(end-p+1:end,:); Y];
X = zeros(T,n*p);
for i=1:p
    X(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end
X = [ones(T,1) X];
k_beta = n*(n-1)/2;         %#ok<NASGU>
k_alp = n^2*p + n;
k = k_alp/n;                %#ok<NASGU>

kappa = .2^2;
kappa1 = .2^2; kappa3 = 100; kappa4 = .2^2;
if is_kappasym
    kappa2 = kappa1;        %#ok<NASGU>
else
    kappa2 = (.2^2)^2;      %#ok<NASGU>
end
sig2_hat = []; U_hat = []; Hyper = struct();  %#ok<NASGU>
switch imodel
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

lml = []; lmlstd = []; store_w = [];
store_Sig = []; store_a = []; store_h = []; store_hpara = []; store_kappa = [];
store_alp = []; store_beta = []; store_l = []; store_A = []; store_F = [];
store_o = []; store_po = [];

scripts = {'VAR_NCP', 'VAR_CSV', 'VAR_ARSV_redu', 'VAR_FSV', 'VAR_ARSVO_redu'};
mls = {'', 'ml_var_csv', 'ml_var_arsv_redu', 'ml_var_fsv', 'ml_var_arsvo_redu'};
resolved = which(scripts{imodel});
assert(strncmpi(resolved, tmp, numel(tmp)), ...
    '%s must resolve from the tempdir copy, got %s', scripts{imodel}, resolved);
resolved = which(mls{imodel});
assert(strncmpi(resolved, tmp, numel(tmp)), ...
    '%s must resolve from the tempdir copy, got %s', mls{imodel}, resolved);

rng(seed, 'twister');
switch imodel
    case 2, VAR_CSV;
    case 3, VAR_ARSV_redu;
    case 4, VAR_FSV;
    case 5, VAR_ARSVO_redu;
end
s = rng;

out = struct('lml',lml, 'lmlstd',lmlstd, 'store_w',store_w, ...
    'store_Sig',store_Sig, 'store_a',store_a, 'store_h',store_h, ...
    'store_hpara',store_hpara, 'store_kappa',store_kappa, 'store_alp',store_alp, ...
    'store_beta',store_beta, 'store_l',store_l, 'store_A',store_A, ...
    'store_F',store_F, 'store_o',store_o, 'store_po',store_po, ...
    'model_name',model_name, 'rngstate',s.State);
end
