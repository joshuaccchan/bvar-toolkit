function test_forecast_iterate_springer
% seeded draw-for-draw equivalence of the functionized springer large-BVAR
% FORECAST pipeline with the legacy workspace scripts of
% chan2020_springer_largebvar, primed exactly as main_forecasting.m primes
% them for one real-time vintage t and run at small nsims from tempdir copies:
%   model 2  forecast_BVAR_Minn.m     -> bvar.priors.minn      + iterate('springer_gauss')
%   model 6  forecast_BVAR_CSV.m      -> bvar.priors.niw('largebvar_nc') +
%                                        bvar.sv.csv_armh      + iterate('springer_csv')
%   model 7  forecast_BVAR_CSV_t.m    -> + bvar.sv.nu_studentt + iterate('springer_csv_t')
%   model 8  forecast_BVAR_CSV_t_MA.m -> + inline psi-MH      + iterate('springer_csv_t_ma')
% Model 2 runs at TWO vintages - one with the latest observation missing
% (is_last_miss true: the extra pre-step simulation draw) and one without -
% models 6/7/8 at the missing-latest vintage. Asserts isequal on tmpyhat0,
% tmpyhat1 and the terminal rng state (same randn/rand/gamrnd/iwishrnd call
% sequence). ZERO PATCHES: none of the springer forecast scripts carries a
% clock-seed line (verified below); the tempdir copies run byte-verbatim with
% nsims/burnin overridden from the harness workspace.
%
% The 20-file vintage xlsread (main_forecasting.m lines 34-53) is slow, so it
% runs ONCE and is cached to a .mat inside the tempdir; the vintage assembly
% itself doubles as the first real-vintage verification of
% bvar.forecast.realtime_loaddata against the legacy loaddata.m copy.
%
% Estimation blocks with no core counterpart yet (the per-model conditional
% draws: iwishrnd Sig/A, lam/sigh2 gamrnd, the rho and psi MH steps) are
% kept VERBATIM inside run_core_* below - this test certifies the forecast
% engine; those blocks belong to the springer family pass.
root = getappdata(0, 'bvar_repo_root');
leg = fullfile(root, 'replications', 'chan2020_springer_largebvar', 'legacy');

tmpdir = tempname; mkdir(tmpdir);
ctmp = onCleanup(@() cleanup_tmp(tmpdir));   %#ok<NASGU>
files = {'forecast_BVAR_Minn.m', 'forecast_BVAR_CSV.m', 'forecast_BVAR_CSV_t.m', ...
    'forecast_BVAR_CSV_t_MA.m', 'loaddata.m', 'prior_Minn.m', 'prior_NC.m', ...
    'SURform2.m', 'sample_h.m', 'sample_nu.m', 'llike_CSV_MA.m'};
for kf = 1:numel(files)
    copyfile(fullfile(leg, files{kf}), fullfile(tmpdir, files{kf}));
    assert(isempty(strfind(fileread(fullfile(tmpdir, files{kf})), 'randn(''seed''')), ...
        'unexpected clock-seed line in %s - the zero-patch premise is broken', files{kf}); %#ok<STREMP>
end
addpath(tmpdir);

    % fixed setup (main_forecasting.m lines 26-56, 68-78)
n = 20; p = 4; k = n*p+1; T0 = 41; T = 208;
tcode = [5 5 5 5 1 5 5 1 5 5 5 5 5 1 1 1 1 1 1 5]';
var_type = [1 1 1 1 1 1 3 2 3 3 3 1 1 4 4 4 4 4 4 4];
psi0 = 0; Vpsi = 1;
lpri_psi = @(x) -.5*(x-psi0)^2/Vpsi -1e10*(x<-.99 || x>.99);
nuub = 50;
nuh0 = 5; Sh0 = .01*(nuh0-1);
rho0 = .9; Vrho = .2^2;

    % load all vintages ONCE (xlsread is deprecated-but-working on R2025b and
    % slow), cache to a .mat in the tempdir
cache = fullfile(tmpdir, 'springer_vintages.mat');
if exist(cache, 'file')
    cc = load(cache); rt_data = cc.rt_data; nonrev_data = cc.nonrev_data;
else
    ws = warning('off', 'all');
    cws = onCleanup(@() warning(ws));
    rt_data.var1 = xlsread(fullfile(leg,'ROUTPUTQvQd.xlsx'),'AI70:HA283');
    rt_data.var2 = xlsread(fullfile(leg,'RCONQvQd.xlsx'),'AI70:HA283');
    rt_data.var3 = xlsread(fullfile(leg,'rinvbfQvQd.xlsx'),'AI70:HA283');
    rt_data.var4 = xlsread(fullfile(leg,'rinvresidQvQd.xlsx'),'AI70:HA283');
    rt_data.var5 = xlsread(fullfile(leg,'RNXQvQd.xlsx'),'AI71:HA283');
    rt_data.var6 = xlsread(fullfile(leg,'npiQvQd.xlsx'),'AI70:HA283');
    rt_data.var7 = xlsread(fullfile(leg,'iptMvMd.xlsx'),'EF542:YI1184');
    rt_data.var8 = xlsread(fullfile(leg,'rucQvMd.xlsx'),'AI209:HA848');
    rt_data.var9 = xlsread(fullfile(leg,'employMvMd.xlsx'),'DG302:XJ944');
    rt_data.var10 = xlsread(fullfile(leg,'hMvMd.xlsx'),'AD206:UG848');
    rt_data.var11 = xlsread(fullfile(leg,'hstartsMvMd.xlsx'),'BW206:VX848');
    rt_data.var12 = xlsread(fullfile(leg,'pconQvQd.xlsx'),'AI70:HA283');
    rt_data.var13 = xlsread(fullfile(leg,'pimpQvQd.xlsx'),'AI70:HA283');
    nonrev_data.var14 = xlsread(fullfile(leg,'FEDFUNDS.xls'),'B129:B770');
    nonrev_data.var15 = xlsread(fullfile(leg,'GS1.xls'),'B144:B785');
    nonrev_data.var16 = xlsread(fullfile(leg,'GS10.xls'),'B144:B785');
    nonrev_data.var17 = xlsread(fullfile(leg,'BAAFFM.xls'),'B129:B770');
    nonrev_data.var18 = xlsread(fullfile(leg,'ISM-MAN_PMI.csv'),'B196:B837');
    nonrev_data.var19 = xlsread(fullfile(leg,'ISM-MAN_NEWORDERS.xls'),'B196:B837');
    nonrev_data.var20 = xlsread(fullfile(leg,'SP500.xlsx'),'B62:B705');
    save(cache, 'rt_data', 'nonrev_data');
    clear cws
end

    % find one vintage WITH a missing latest observation and one WITHOUT,
    % verifying core realtime_loaddata against the legacy loaddata copy as we go
t_miss = []; t_nomiss = [];
for t = T0:T
    [dtc, dkc] = bvar.forecast.realtime_loaddata(rt_data, nonrev_data, t, T0, tcode, var_type);
    [dtl, dkl] = loaddata(rt_data, nonrev_data, t, T0, tcode, var_type);   % tempdir legacy copy
    assert(isequaln(dtc, dtl) && isequaln(dkc, dkl), ...
        'realtime_loaddata differs from legacy loaddata at t=%d', t);
    ilm = (sum(isnan(dtc(end, :)), 2) > 0);
    if ilm && isempty(t_miss), t_miss = t; end
    if ~ilm && isempty(t_nomiss), t_nomiss = t; end
    if ~isempty(t_miss) && ~isempty(t_nomiss), break; end
end
assert(~isempty(t_miss) && ~isempty(t_nomiss), ...
    'need both a missing-latest and a complete vintage in %d..%d', T0, t);

D0 = struct('n',n,'p',p,'k',k,'T',T,'T0',T0,'tcode',tcode,'var_type',var_type, ...
    'nuh0',nuh0,'Sh0',Sh0,'rho0',rho0,'Vrho',Vrho,'nuub',nuub,'lpri_psi',lpri_psi);
nsims = 25; burnin = 8; seed = 20260901;

    % {script, core runner, c-hyperparameters, vintage}
cases = { ...
    'forecast_BVAR_Minn',     @run_core_minn,     [.2^2 .1^2 100], t_miss;
    'forecast_BVAR_Minn',     @run_core_minn,     [.2^2 .1^2 100], t_nomiss;
    'forecast_BVAR_CSV',      @run_core_csv,      [.2^2 100 NaN],  t_miss;
    'forecast_BVAR_CSV_t',    @run_core_csv_t,    [.2^2 100 NaN],  t_miss;
    'forecast_BVAR_CSV_t_MA', @run_core_csv_t_ma, [.2^2 100 NaN],  t_miss};
for kc = 1:size(cases, 1)
    [script, corefun, cvec, t] = cases{kc, :};
    D = D0;
    D.c1 = cvec(1); D.c2 = cvec(2); D.c3 = cvec(3);
    D = add_vintage(D, rt_data, nonrev_data, t);
    L = run_legacy(script, tmpdir, D, nsims, burnin, seed);
    R = corefun(D, nsims, burnin, seed);
    lbl = sprintf('%s @ t=%d (is_last_miss=%d)', script, t, D.is_last_miss);
    assert(isequal(L.tmpyhat0, R.tmpyhat0), '%s: tmpyhat0 differs', lbl);
    assert(isequal(L.tmpyhat1, R.tmpyhat1), '%s: tmpyhat1 differs', lbl);
    assert(isequal(L.rngstate, R.rngstate), '%s: rng call sequence differs', lbl);
end
end

% ---------------------------------------------------------------------------
function cleanup_tmp(tmpdir)
if any(strcmpi(strsplit(path, pathsep), tmpdir))
    rmpath(tmpdir);
end
if exist(tmpdir, 'dir')
    rmdir(tmpdir, 's');
end
end

function D = add_vintage(D, rt_data, nonrev_data, t)
% per-vintage assembly, verbatim main_forecasting.m lines 104-119 (both sides
% receive the SAME precomputed panel; the assembly itself has no rng calls)
[data_t, data_tpk] = bvar.forecast.realtime_loaddata(rt_data, nonrev_data, ...
    t, D.T0, D.tcode, D.var_type);
is_last_miss = (sum(isnan(data_t(end, :)), 2) > 0);
if is_last_miss
    data_t = data_t(1:end-1, :);
end
ind = sum(isnan(data_t), 2);
ind_last = find(ind > 0, 1, 'last');
if ~isempty(ind_last)
    data_t = data_t(ind_last+1:end, :);
end
D.Y0 = data_t(1:4, :);
D.shortYt = data_t(5:end, :);
D.Tt = size(D.shortYt, 1);
D.data_tpk = data_tpk;
D.is_last_miss = is_last_miss;
D.t = t;
end

function out = run_legacy(script, tmpdir, D, nsims, burnin, seed) %#ok<INUSD> % nsims/burnin are read by the dispatched script from this workspace
% unpack EVERY workspace variable the legacy forecast scripts read, with the
% small nsims/burnin already assigned, then dispatch the byte-verbatim script
n = D.n; p = D.p; k = D.k; T = D.T; t = D.t;                    %#ok<NASGU>
Y0 = D.Y0; shortYt = D.shortYt; Tt = D.Tt;                      %#ok<NASGU>
data_tpk = D.data_tpk; is_last_miss = D.is_last_miss;           %#ok<NASGU>
c1 = D.c1; c2 = D.c2; c3 = D.c3;                                %#ok<NASGU>
nuh0 = D.nuh0; Sh0 = D.Sh0; rho0 = D.rho0; Vrho = D.Vrho;       %#ok<NASGU>
nuub = D.nuub; lpri_psi = D.lpri_psi;                           %#ok<NASGU>

tmpyhat0 = []; tmpyhat1 = [];   % pre-declare script outputs (parser binding)

resolved = which(script);
assert(strncmpi(resolved, tmpdir, numel(tmpdir)), ...
    '%s must resolve from the tempdir copy, got %s', script, resolved);

rng(seed, 'twister');
switch script
    case 'forecast_BVAR_Minn',      forecast_BVAR_Minn;
    case 'forecast_BVAR_CSV',       forecast_BVAR_CSV;
    case 'forecast_BVAR_CSV_t',     forecast_BVAR_CSV_t;
    case 'forecast_BVAR_CSV_t_MA',  forecast_BVAR_CSV_t_MA;
    otherwise, error('unknown legacy script %s', script);
end
s = rng;
out = struct('tmpyhat0', tmpyhat0, 'tmpyhat1', tmpyhat1, 'rngstate', s.State);
end

function cfg = fc_cfg(D)
cfg = struct('shortYt', D.shortYt, 'data_tpk', D.data_tpk, ...
    'is_last_miss', D.is_last_miss, 'p', D.p, 't', D.t, 'T', D.T);
end

function out = run_core_minn(D, nsims, burnin, seed)
% functionized model 2 [forecast_BVAR_Minn.m]: bvar.priors.minn +
% bvar.util.surform2, forecast via iterate('springer_gauss')
n = D.n; p = D.p; k = D.k; Y0 = D.Y0; shortYt = D.shortYt; Tt = D.Tt;
rng(seed, 'twister');
Yt = reshape(shortYt', n*Tt, 1);
tmpyhat0 = zeros(nsims, 2*n+1);
tmpyhat1 = zeros(nsims, 2*n+1);
[beta_Minn, V_Minn, Sig_hat] = bvar.priors.minn(p, D.c1, D.c2, D.c3, Y0, shortYt, p);  % legacy prior_Minn
[~, Z] = bvar.util.build_lags([Y0(end-p+1:end, :); shortYt], p);
X = bvar.util.surform2(Z, n);                            % legacy SURform2
XiSig = X'*kron(speye(Tt), sparse(1:n, 1:n, 1./Sig_hat));
Kbeta = sparse(1:n*k, 1:n*k, 1./V_Minn) + XiSig*X;
C_Kbeta = chol(Kbeta, 'lower');
beta_hat = C_Kbeta'\(C_Kbeta\(beta_Minn./V_Minn + XiSig * Yt));
CSig = sparse(1:n, 1:n, sqrt(Sig_hat));                 % SPARSE, as in the legacy script
cfg = fc_cfg(D);
for isim = 1:nsims + burnin
    beta = beta_hat + C_Kbeta'\randn(k*n, 1);
    if isim > burnin
        isave = isim - burnin;
        A = reshape(beta, k, n);                        % legacy line 35 (caller-side)
        fcr = bvar.forecast.iterate('springer_gauss', ...
            struct('A', A, 'CSig', CSig, 'dSig', Sig_hat'), cfg);
        tmpyhat0(isave, :) = fcr(1, :);
        tmpyhat1(isave, :) = fcr(2, :);
    end
end
s = rng;
out = struct('tmpyhat0', tmpyhat0, 'tmpyhat1', tmpyhat1, 'rngstate', s.State);
end

function out = run_core_csv(D, nsims, burnin, seed)
% functionized model 6 [forecast_BVAR_CSV.m]: bvar.priors.niw('largebvar_nc') +
% bvar.sv.csv_armh; the Sig/A, sigh2 and rho-MH draws stay verbatim inline
% (springer family pass); forecast via iterate('springer_csv')
n = D.n; p = D.p; k = D.k; Y0 = D.Y0; shortYt = D.shortYt; Tt = D.Tt;
nuh0 = D.nuh0; Sh0 = D.Sh0; rho0 = D.rho0; Vrho = D.Vrho;
rng(seed, 'twister');
tmpyhat0 = zeros(nsims, 2*n+1);
tmpyhat1 = zeros(nsims, 2*n+1);
[A0, VA0, nu0, S0] = bvar.priors.niw(p, [D.c1 D.c2], Y0, shortYt, 'largebvar_nc');  % legacy prior_NC
[~, Z] = bvar.util.build_lags([Y0(end-p+1:end, :); shortYt], p);
h = zeros(Tt, 1);
rho = .8;
sigh2 = .1;
cfg = fc_cfg(D);
for isim = 1:nsims + burnin
        % sample Sig and A [verbatim lines 31-40]
    iOh = sparse(1:Tt, 1:Tt, exp(-h));
    ZiOh = Z'*iOh;
    KA = sparse(1:k, 1:k, 1./VA0) + ZiOh*Z;
    Ahat = KA\(sparse(1:k, 1:k, VA0)\A0 + ZiOh*shortYt);
    Shat = S0 + A0'*sparse(1:k, 1:k, 1./VA0)*A0 + shortYt'*iOh*shortYt ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2;
    Sig = iwishrnd(Shat, nu0+Tt);
    CSig = chol(Sig, 'lower');
    A = Ahat + (chol(KA, 'lower')'\randn(k, n))*CSig';
        % sample h [lines 43-46] via the extracted core block
    U = shortYt - Z*A;
    tmph = (U/CSig');
    s2_h = sum(tmph.^2, 2);
    h = bvar.sv.csv_armh(s2_h, rho, sigh2, h, n);        % legacy sample_h
        % sample sigh2 [lines 49-50]
    eh = [h(1)*sqrt(1-rho^2);  h(2:end)-rho*h(1:end-1)];
    sigh2 = 1/gamrnd(nuh0+Tt/2, 1/(Sh0 + sum(eh.^2)/2));
        % sample rho [lines 53-63; CSV bound .999]
    Krho = 1/Vrho + sum(h(1:Tt-1).^2)/sigh2;
    rhohat = Krho\(rho0/Vrho + h(1:Tt-1)'*h(2:Tt)/sigh2);
    rhoc = rhohat + sqrt(Krho)'\randn;
    grho = @(x) -.5*log(sigh2./(1-x.^2))-.5*(1-x.^2)/sigh2*h(1)^2;
    if abs(rhoc) < .999
        alpMH = exp(grho(rhoc)-grho(rho));
        if alpMH > rand
            rho = rhoc;
        end
    end
    if isim > burnin
        isave = isim - burnin;
        fcr = bvar.forecast.iterate('springer_csv', ...
            struct('A', A, 'CSig', CSig, 'Sig', Sig, 'h', h, 'rho', rho, 'sigh2', sigh2), cfg);
        tmpyhat0(isave, :) = fcr(1, :);
        tmpyhat1(isave, :) = fcr(2, :);
    end
end
s = rng;
out = struct('tmpyhat0', tmpyhat0, 'tmpyhat1', tmpyhat1, 'rngstate', s.State);
end

function out = run_core_csv_t(D, nsims, burnin, seed)
% functionized model 7 [forecast_BVAR_CSV_t.m]: adds the lam gamrnd draw and
% bvar.sv.nu_studentt (legacy sample_nu); rho bound .99 here (CSV uses .999);
% forecast via iterate('springer_csv_t')
n = D.n; p = D.p; k = D.k; Y0 = D.Y0; shortYt = D.shortYt; Tt = D.Tt;
nuh0 = D.nuh0; Sh0 = D.Sh0; rho0 = D.rho0; Vrho = D.Vrho; nuub = D.nuub;
rng(seed, 'twister');
tmpyhat0 = zeros(nsims, 2*n+1);
tmpyhat1 = zeros(nsims, 2*n+1);
[A0, VA0, nu0, S0] = bvar.priors.niw(p, [D.c1 D.c2], Y0, shortYt, 'largebvar_nc');
[~, Z] = bvar.util.build_lags([Y0(end-p+1:end, :); shortYt], p);
h = zeros(Tt, 1);
nu = 5;
rho = .8;
sigh2 = .1;
lam = 1./gamrnd(nu/2, 2/nu, Tt, 1);
cfg = fc_cfg(D);
for isim = 1:nsims + burnin
        % sample Sig and A [verbatim lines 31-41]
    iOm = sparse(1:Tt, 1:Tt, exp(-h)./lam);
    ZiOm = Z'*iOm;
    KA = sparse(1:k, 1:k, 1./VA0) + ZiOm*Z;
    Ahat = KA\(sparse(1:k, 1:k, VA0)\A0 + ZiOm*shortYt);
    Shat = S0 + A0'*sparse(1:k, 1:k, 1./VA0)*A0 + shortYt'*iOm*shortYt ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2;
    Sig = iwishrnd(Shat, nu0+Tt);
    CSig = chol(Sig, 'lower');
    A = Ahat + (chol(KA, 'lower')'\randn(k, n))*CSig';
        % sample h [lines 44-47]
    U = shortYt - Z*A;
    tmph = (U/CSig');
    s2_h = sum(tmph.^2, 2)./lam;
    h = bvar.sv.csv_armh(s2_h, rho, sigh2, h, n);
        % sample lam [lines 50-53]
    U = shortYt - Z*A;
    tmph = (U/CSig');
    s2 = sum(tmph.^2, 2)./exp(h);
    lam = 1./gamrnd((n+nu)/2, 2./(s2+nu));
        % sample nu [line 56]
    nu = bvar.sv.nu_studentt(lam, nu, nuub);             % legacy sample_nu
        % sample sigh2 [lines 59-60]
    eh = [h(1)*sqrt(1-rho^2);  h(2:end)-rho*h(1:end-1)];
    sigh2 = 1/gamrnd(nuh0+Tt/2, 1/(Sh0 + sum(eh.^2)/2));
        % sample rho [lines 63-73; bound .99]
    Krho = 1/Vrho + sum(h(1:Tt-1).^2)/sigh2;
    rhohat = Krho\(rho0/Vrho + h(1:Tt-1)'*h(2:Tt)/sigh2);
    rhoc = rhohat + sqrt(Krho)'\randn;
    grho = @(x) -.5*log(sigh2./(1-x.^2))-.5*(1-x.^2)/sigh2*h(1)^2;
    if abs(rhoc) < .99
        alpMH = exp(grho(rhoc)-grho(rho));
        if alpMH > rand
            rho = rhoc;
        end
    end
    if isim > burnin
        isave = isim - burnin;
        fcr = bvar.forecast.iterate('springer_csv_t', ...
            struct('A', A, 'CSig', CSig, 'Sig', Sig, 'h', h, 'rho', rho, ...
            'sigh2', sigh2, 'nu', nu), cfg);
        tmpyhat0(isave, :) = fcr(1, :);
        tmpyhat1(isave, :) = fcr(2, :);
    end
end
s = rng;
out = struct('tmpyhat0', tmpyhat0, 'tmpyhat1', tmpyhat1, 'rngstate', s.State);
end

function out = run_core_csv_t_ma(D, nsims, burnin, seed)
% functionized model 8 [forecast_BVAR_CSV_t_MA.m]: adds the MA(1) psi-MH
% (verbatim inline, incl. fminunc/fminbnd mode-finding; llike_CSV_MA resolves
% to the tempdir LEGACY copy - the springer/realtime copies omit the
% -n/2*sum(h) term, see the never-merge list); forecast via
% iterate('springer_csv_t_ma')
n = D.n; p = D.p; k = D.k; Y0 = D.Y0; shortYt = D.shortYt; Tt = D.Tt;
nuh0 = D.nuh0; Sh0 = D.Sh0; rho0 = D.rho0; Vrho = D.Vrho; nuub = D.nuub;
lpri_psi = D.lpri_psi;
rng(seed, 'twister');
tmpyhat0 = zeros(nsims, 2*n+1);
tmpyhat1 = zeros(nsims, 2*n+1);
[A0, VA0, nu0, S0] = bvar.priors.niw(p, [D.c1 D.c2], Y0, shortYt, 'largebvar_nc');
[~, Z] = bvar.util.build_lags([Y0(end-p+1:end, :); shortYt], p);
h = zeros(Tt, 1);
nu = 5;
rho = .8;
sigh2 = .1;
lam = 1./gamrnd(nu/2, 2/nu, Tt, 1);
psi = .1;
Hpsi = speye(Tt) + psi*sparse(2:Tt, 1:(Tt-1), ones(1, Tt-1), Tt, Tt);
options = optimset('Display', 'off', 'LargeScale', 'off');
psihat = psi;
cfg = fc_cfg(D);
cfg.Z = Z;
for isim = 1:nsims + burnin
        % sample Sig and A [verbatim lines 34-46]
    Ztld = Hpsi\Z;
    Ytld = Hpsi\shortYt;
    iO_h_lam_psi = sparse(1:Tt, 1:Tt, [1/(1+psi^2)*exp(-h(1)); exp(-h(2:end))]./lam);
    ZiO = Ztld'*iO_h_lam_psi;
    KA = sparse(1:k, 1:k, 1./VA0) + ZiO*Ztld;
    Ahat = KA\(sparse(1:k, 1:k, VA0)\A0 + ZiO*Ytld);
    Shat = S0 + A0'*sparse(1:k, 1:k, 1./VA0)*A0 + Ytld'*iO_h_lam_psi*Ytld ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2;
    Sig = iwishrnd(Shat, nu0+Tt);
    CSig = chol(Sig, 'lower');
    A = Ahat + (chol(KA, 'lower')'\randn(k, n))*CSig';
        % sample lam [lines 49-53]
    U = shortYt - Z*A;
    Utld = Hpsi\U;
    tmph = (Utld/CSig');
    s2_lam = sum(tmph.^2, 2)./exp(h);
    lam = 1./gamrnd((n+nu)/2, 2./(s2_lam+nu));
        % sample h [lines 56-58]
    s2_h = sum(tmph.^2, 2)./lam;
    s2_h(1) = s2_h(1)/(1+psi^2);
    h = bvar.sv.csv_armh(s2_h, rho, sigh2, h, n);
        % sample sigh2 [lines 61-62]
    eh = [h(1)*sqrt(1-rho^2);  h(2:end)-rho*h(1:end-1)];
    sigh2 = 1/gamrnd(nuh0+Tt/2, 1/(Sh0 + sum(eh.^2)/2));
        % sample rho [lines 65-75; bound .99]
    Krho = 1/Vrho + sum(h(1:Tt-1).^2)/sigh2;
    rhohat = Krho\(rho0/Vrho + h(1:Tt-1)'*h(2:Tt)/sigh2);
    rhoc = rhohat + sqrt(Krho)'\randn;
    grho = @(x) -.5*log(sigh2./(1-x.^2))-.5*(1-x.^2)/sigh2*h(1)^2;
    if abs(rhoc) < .99
        alpMH = exp(grho(rhoc)-grho(rho));
        if alpMH > rand
            rho = rhoc;
        end
    end
        % sample psi [verbatim lines 78-102]
    U_psi = U./repmat(sqrt(lam), 1, n);
    lp_psi = @(x) llike_CSV_MA(x, U_psi, Sig, h) + lpri_psi(x);
    if (mod(isim, 100) == 0) || isim == 1
        [psihat, fval, exitflag, output, grad, hess] ...
            = fminunc(@(x)-lp_psi(x), psihat, options);  %#ok<ASGLU>
        [tmpCpsi, flag] = chol(hess, 'lower');           %#ok<ASGLU>
        if flag == 0
            Kpsic = hess;
        else
            Kpsic = 1/.05^2;
        end
    else
        psihat = fminbnd(@(x)-lp_psi(x), -.99, .99);
    end
    psic = psihat + 1/sqrt(Kpsic)*randn;
    if abs(psic) < .99
        alpMH =  lp_psi(psic) - lp_psi(psi) + ...
            -.5*(psi-psihat)^2*Kpsic + .5*(psic-psihat)^2*Kpsic;
    else
        alpMH = -inf;
    end
    if alpMH > log(rand)
        psi = psic;
        Hpsi = speye(Tt) + psi*sparse(2:Tt, 1:(Tt-1), ones(1, Tt-1), Tt, Tt);
    end
        % sample nu [line 105]
    nu = bvar.sv.nu_studentt(lam, nu, nuub);
    if isim > burnin
        isave = isim - burnin;
        fcr = bvar.forecast.iterate('springer_csv_t_ma', ...
            struct('A', A, 'CSig', CSig, 'Sig', Sig, 'h', h, 'rho', rho, ...
            'sigh2', sigh2, 'nu', nu, 'psi', psi, 'Hpsi', Hpsi), cfg);
        tmpyhat0(isave, :) = fcr(1, :);
        tmpyhat1(isave, :) = fcr(2, :);
    end
end
s = rng;
out = struct('tmpyhat0', tmpyhat0, 'tmpyhat1', tmpyhat1, 'rngstate', s.State);
end
