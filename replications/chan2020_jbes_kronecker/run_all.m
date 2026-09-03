% chan2020_jbes_kronecker/run_all - functionized full-sample estimation
% pipeline of Chan (2020, JBES): the eight large BVARs with Kronecker error
% covariance structures dispatched by legacy main_BVAR.m -
%   1 'BVAR'         (analytic NIW posterior, no MCMC)      BVAR.m
%   2 'BVAR-t'       (multivariate-t errors)                BVAR_t.m
%   3 'BVAR-CSV'     (common stochastic volatility)         BVAR_CSV.m
%   4 'BVAR-MA'      (MA(1) errors)                         BVAR_MA.m
%   5 'BVAR-t-CSV'   (t + CSV)                              BVAR_t_CSV.m
%   6 'BVAR-t-MA'    (t + MA(1))                            BVAR_t_MA.m
%   7 'BVAR-CSV-MA'  (CSV + MA(1))                          BVAR_CSV_MA.m
%   8 'BVAR-CSV-t-MA' (CSV + t + MA(1))                     BVAR_CSV_t_MA.m
%
%   out = run_all(model, nsim, burnin, seed)
%
%   model  - 1..8 or one of the names above (default preset model_default = 8,
%            the legacy main_BVAR.m line 19 setting)
%   nsim   - posterior draws kept (default preset nsims_default = 30000)
%   burnin - burn-in sweeps (default preset burnin_default = 5000)
%   seed   - optional; when nonempty, rng(seed,'twister') is set FIRST, before
%            any draw. When omitted/empty the ambient rng state is used as-is.
%
% Functionized 2026-09-02 (step 8, Kronecker family pass, part 1: estimation +
% marginal likelihood; the realtime forecasting drivers are part 2). Marginal
% likelihoods are NOT computed here - run_ml.m in this folder runs this
% estimation and then the extracted bvt.ml.kron_bvar_* computation on its
% output (exposing the bugcompat option for the two ml scripts with known
% legacy defects; see run_ml and tests/variant_map.md).
%
% Reproduces the legacy pipeline main_BVAR.m -> BVAR*.m draw-for-draw bitwise
% (verified by tests/unit/test_kron_equivalence.m at small nsim) with ONE
% deliberate divergence: every legacy MCMC script (models 2-8) re-seeds the
% global stream from the wall clock at its "MCMC starts here" banner
% (randn('seed',sum(clock*100)); rand('seed',sum(clock*1000))), which makes
% as-shipped runs irreproducible AND switches MATLAB to the legacy v4/v5
% generators; run_all drops that line so the caller controls seeding via
% `seed` (or the ambient state) on the modern twister stream - same rationale
% and mechanics as the step-5 MAHP and step-7 OISV notes in
% tests/variant_map.md. Model 2/5/6/8 chain-init gamrnd draws sit BEFORE the
% legacy clock-seed line, so under the patched semantics the whole run
% (chain init included) consumes one coherent stream in script order. The
% legacy wall-clock timing displays and all figure windows (imagesc heat
% maps, density plots, histograms) are not reproduced; the disp banners and
% the mod-5000 loop counter are.
%
% All constants come from preset.m in this folder (each field cites its legacy
% source line); the data file is read from legacy/ READ-ONLY; core reuse:
%   bvt.priors.niw('kron_script')  = legacy construct_prior_A (+ the S0/nu0
%                                    assignments its callers make first),
%   bvt.util.build_lags            = the inline X construction,
%   bvt.sv.csv_armh                = legacy sample_h (root copy),
%   bvt.sv.nu_studentt             = legacy sample_nu (root copy; normpdf-form
%                                    MH ratio - see that header's note),
%   bvt.ml.llike_ma / llike_csv_ma = legacy llike_MA / ROOT llike_CSV_MA
%                                    inside the psi-MH steps.
% The Sig/A joint draw, lam, sigh2, rho-MH and psi-MH blocks have no core
% counterpart and live verbatim in the per-model subfunctions below.
%
% Output struct: the legacy stores and script-tail posterior summaries (per
% model), acceptance counters (legacy names), the design (shortY, X, Y0,
% T/n/p/k), the prior struct `pri` (A0, VA0, nu0, S0 + the model's extra
% hyperparameters), the preset used, and `state` = the final workspace values
% the legacy ml_* scripts consume as leftovers (final chain draws of
% Sig/A/h/lam/rho/sigh2/psi/nu, the final psi-MH mode psihat, counters) -
% made explicit so run_ml/bvt.ml.* can reproduce the legacy ML computations
% exactly (including their leftover-workspace defects under bugcompat).
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics,
% 38(1), 68-79.

function out = run_all(model, nsim, burnin, seed)
thisdir = fileparts(mfilename('fullpath'));

    % make bvt.* resolvable when called standalone
if isempty(which('bvt.ml.lniwpdf'))
    root = fileparts(fileparts(thisdir));
    addpath(fullfile(root, 'core'));
end

    % constants: preset.m in THIS folder (cd guard pins name resolution)
od = cd(thisdir);
guard = onCleanup(@() cd(od));
pr = preset();
clear guard

if nargin < 1 || isempty(model),  model  = pr.model_default;   end
if nargin < 2 || isempty(nsim),   nsim   = pr.nsims_default;   end
if nargin < 3 || isempty(burnin), burnin = pr.burnin_default;  end
if nargin < 4, seed = []; end
names = {'BVAR','BVAR-t','BVAR-CSV','BVAR-MA','BVAR-t-CSV','BVAR-t-MA', ...
    'BVAR-CSV-MA','BVAR-CSV-t-MA'};
if ischar(model) || isstring(model)
    im = find(strcmpi(model, names), 1);
    if isempty(im)
        error('run_all:model', 'unknown model ''%s''', model);
    end
    model = im;
elseif ~(isscalar(model) && any(model == 1:8))
    error('run_all:model', 'model must be 1..8 or a model name');
end
if ~isempty(seed)
    rng(seed, 'twister');
end

    % data (legacy folder, read-only) and design [main_BVAR.m lines 26-31]
p = pr.p;
data_Q = load(fullfile(thisdir, 'legacy', pr.data_file));
data = data_Q(:, pr.var_cols);
Y0 = data(1:pr.n0, :);                  % save the first 4 obs as the initial conditions
shortY = data(pr.n0+1:end, :);
[T, n] = size(shortY);
k = n*p+1;                              % # of coefficients in each equation

    % prior [each BVAR*.m prior block: S0 = eye(n); nu0 = n+3;
    % construct_prior_A -> bvt.priors.niw 'kron_script' variant]
[A0, VA0, nu0, S0] = bvt.priors.niw(p, pr.minn_kappa, Y0, shortY, 'kron_script');
pri = struct('A0', A0, 'VA0', VA0, 'nu0', nu0, 'S0', S0);

    % X [each BVAR*.m "construct X" block; identical to the legacy inline loop]
[~, X] = bvt.util.build_lags([Y0(end-p+1:end, :); shortY], p);

switch model
    case 1, res = post_bvar(shortY, X, pri, T, n, k);
    case 2, res = mcmc_t(shortY, X, pri, T, n, k, nsim, burnin, pr);
    case 3, res = mcmc_csv(shortY, X, pri, T, n, k, nsim, burnin, pr);
    case 4, res = mcmc_ma(shortY, X, pri, T, n, k, nsim, burnin, pr);
    case 5, res = mcmc_t_csv(shortY, X, pri, T, n, k, nsim, burnin, pr);
    case 6, res = mcmc_t_ma(shortY, X, pri, T, n, k, nsim, burnin, pr);
    case 7, res = mcmc_csv_ma(shortY, X, pri, T, n, k, nsim, burnin, pr);
    case 8, res = mcmc_csv_t_ma(shortY, X, pri, T, n, k, nsim, burnin, pr);
end

out = res;
out.model = model;
out.model_name = names{model};
out.nsims = nsim;
out.burnin = burnin;
out.seed = seed;
out.T = T; out.n = n; out.p = p; out.k = k;
out.shortY = shortY; out.X = X; out.Y0 = Y0;
out.pri = res.pri;
out.preset = pr;
end

% -------------------------------------------------------------------------
function res = post_bvar(shortY, X, pri, T, ~, k)
% BVAR.m functionized (lines 13-26; the figure block and the inline cp_ml
% marginal-likelihood block are not reproduced here - the latter is
% bvt.ml.kron_bvar, called by run_ml). No MCMC, no rng.
A0 = pri.A0; VA0 = pri.VA0; nu0 = pri.nu0; S0 = pri.S0;

XX = X'*X;
Atilde = XX\(X'*shortY);
KA = sparse(1:k,1:k,1./VA0) + XX;
    % posterior mean of the VAR coefficients, arranged as a k by n matrix
Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XX*Atilde);
Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + shortY'*shortY ...
    - Ahat'*KA*Ahat;
Shat = (Shat+Shat')/2;

res = struct('Atilde', Atilde, 'KA', KA, 'Ahat', Ahat, 'Shat', Shat);
res.Sig_hat = Shat/(T+nu0);             % the (A,Sig) evaluation pair of the ML block
res.pri = pri;
res.state = struct();
end

% -------------------------------------------------------------------------
function res = mcmc_t(shortY, X, pri, T, n, k, nsims, burnin, pr)
% BVAR_t.m functionized line-for-line (clock-seed line 34 dropped; figures
% not reproduced). Draw order per sweep: (Sig,A) -> lam -> nu.
A0 = pri.A0; VA0 = pri.VA0; nu0 = pri.nu0; S0 = pri.S0;
nuub = pr.nuub;                         % BVAR_t.m line 11
pri.nuub = nuub;

    % initialize for storage [lines 21-26]
store_Sig = zeros(n,n);
store_A = zeros(k,n);
store_nu = zeros(nsims,1);
store_lam = zeros(nsims,T);
nugrid = linspace(2,nuub,700)';
store_pnu = zeros(700,1);

    % initialize the chain [lines 29-31]
nu = 5;
lam = 1./gamrnd(nu/2,2/nu,T,1);
countnu = 0;

    % (legacy clock-seed line 34 deliberately dropped)
disp('Starting MCMC for BVAR-t.... ');

for isim = 1:nsims + burnin

        %% sample Sig and A
    iOm = sparse(1:T,1:T,1./lam);
    XiOm = X'*iOm;
    KA = sparse(1:k,1:k,1./VA0) + XiOm*X;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiOm*shortY);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + shortY'*iOm*shortY ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    Sig = iwishrnd(Shat,nu0+T);
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig';

        %% sample lam
    U = shortY - X*A;
    tmp = (U/CSig');
    s2 = sum(tmp.^2,2);
    lam = 1./gamrnd((n+nu)/2,2./(s2+nu));

        %% sample nu [legacy sample_nu -> bvt.sv.nu_studentt]
    [nu,flag,fnu] = bvt.sv.nu_studentt(lam,nu,nuub);
    countnu = countnu + flag;

    if isim > burnin
        isave = isim - burnin;
        store_A = store_A + A;
        store_Sig = store_Sig + Sig;
        store_nu(isave,:) = nu;
        store_lam(isave,:) = lam';

        % compute the density of nu
        tmpden = fnu(nugrid);
        tmpden = exp(tmpden-max(tmpden));
        tmpden = tmpden/(sum(tmpden)*(nugrid(2)-nugrid(1)));
        store_pnu = store_pnu + tmpden;
    end

    if ( mod(isim, pr.progress_every) ==0 )
        disp(  [ num2str(isim) ' loops... ' ] )
    end

end

    % script-tail posterior summaries [lines 85-88]
res = struct();
res.A_mean = store_A/nsims;
res.Sig_mean = store_Sig/nsims;
res.nu_mean = mean(store_nu)';
res.pnu_mean = store_pnu/nsims;
res.store_A = store_A; res.store_Sig = store_Sig;
res.store_nu = store_nu; res.store_lam = store_lam;
res.store_pnu = store_pnu;
res.countnu = countnu;
res.pri = pri;
res.state = struct('nu', nu, 'lam', lam, 'Sig', Sig, 'A', A, 'countnu', countnu);
end

% -------------------------------------------------------------------------
function res = mcmc_csv(shortY, X, pri, T, n, k, nsims, burnin, pr)
% BVAR_CSV.m functionized line-for-line (clock-seed line 35 dropped; figures
% not reproduced). Draw order per sweep: (Sig,A) -> h -> sigh2 -> rho.
A0 = pri.A0; VA0 = pri.VA0; nu0 = pri.nu0; S0 = pri.S0;
nuh0 = pr.nuh0; Sh0 = pr.Sh0;           % BVAR_CSV.m line 11
rho0 = pr.rho0; Vrho = pr.Vrho;         % line 12
pri.nuh0 = nuh0; pri.Sh0 = Sh0; pri.rho0 = rho0; pri.Vrho = Vrho;

    % initialize for storage [lines 22-25]
store_Sig = zeros(n,n);
store_A = zeros(k,n);
store_h = zeros(nsims,T);
store_theta = zeros(nsims,2);

    % initialize the chain [lines 28-32]
h = zeros(T,1);
rho = .8;
sigh2 = .1;
Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T); %#ok<NASGU> % dead in the loop (sample_h rebuilds it); kept verbatim
counth = 0; countrho = 0;

    % (legacy clock-seed line 35 deliberately dropped)
disp('Starting MCMC for BVAR-CSV.... ');

for isim = 1:nsims + burnin

        %% sample Sig and A
    iOh = sparse(1:T,1:T,exp(-h));
    XiOh = X'*iOh;
    KA = sparse(1:k,1:k,1./VA0) + XiOh*X;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiOh*shortY);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + shortY'*iOh*shortY ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    Sig = iwishrnd(Shat,nu0+T);
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig';

        %% sample h [legacy sample_h -> bvt.sv.csv_armh]
    U = shortY - X*A;
    tmp = (U/CSig');
    s2 = sum(tmp.^2,2);
    [h, flag] = bvt.sv.csv_armh(s2,rho,sigh2,h,n);
    counth = counth + flag;

        %% sample sigh2
    eh = [h(1)*sqrt(1-rho^2);  h(2:end)-rho*h(1:end-1)];
    sigh2 = 1/gamrnd(nuh0+T/2,1/(Sh0 + sum(eh.^2)/2));

        %% sample rho [estimation truncation bound .9999 - see preset]
    Krho = 1/Vrho + sum(h(1:T-1).^2)/sigh2;
    rhohat = Krho\(rho0/Vrho + h(1:T-1)'*h(2:T)/sigh2);
    rhoc = rhohat + sqrt(Krho)'\randn;
    grho = @(x) -.5*log(sigh2./(1-x.^2))-.5*(1-x.^2)/sigh2*h(1)^2;
    if abs(rhoc)<.9999
        alpMH = exp(grho(rhoc)-grho(rho));
        if alpMH>rand
            rho = rhoc;
            countrho = countrho+1;
            Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T); %#ok<NASGU>
        end
    end

    if isim > burnin
        isave = isim - burnin;
        store_A = store_A + A;
        store_Sig = store_Sig + Sig;
        store_h(isave,:) = h';
        store_theta(isave,:) = [rho sigh2];
    end

    if ( mod(isim, pr.progress_every) ==0 )
        disp(  [ num2str(isim) ' loops... ' ] )
    end

end

    % script-tail posterior summaries [lines 95-99]
res = struct();
res.A_mean = store_A/nsims;
res.Sig_mean = store_Sig/nsims;
res.h_mean = mean(store_h)';
res.theta_mean = mean(store_theta)';
res.CSV_std_mean = mean(exp(store_h/2))';
res.store_A = store_A; res.store_Sig = store_Sig;
res.store_h = store_h; res.store_theta = store_theta;
res.counth = counth; res.countrho = countrho;
res.pri = pri;
res.state = struct('h', h, 'rho', rho, 'sigh2', sigh2, 'Sig', Sig, 'A', A, ...
    'counth', counth, 'countrho', countrho);
end

% -------------------------------------------------------------------------
function res = mcmc_ma(shortY, X, pri, T, n, k, nsims, burnin, pr)
% BVAR_MA.m functionized line-for-line (clock-seed line 38 dropped; figures
% not reproduced). Draw order per sweep: (Sig,A) -> psi.
A0 = pri.A0; VA0 = pri.VA0; nu0 = pri.nu0; S0 = pri.S0;
psi0 = pr.psi0; Vpsi = pr.Vpsi;         % BVAR_MA.m line 9
lpri_psi = @(x) -.5*(x-psi0)^2/Vpsi -1e10*(x<-.99 || x>.99);    % line 10
pri.psi0 = psi0; pri.Vpsi = Vpsi;

    % initialize for storage [lines 22-28]
store_Sig = zeros(n,n);
store_A = zeros(k,n);
store_psi = zeros(nsims,1);
ngrid = 300;
psigrid = linspace(-.1,.5,ngrid)';      %#ok<NASGU> % legacy line 26; read only by the commented-out ppsi block
store_ppsi = zeros(ngrid,1);
countpsi = 0;

    % initialize the chain [lines 31-35]
psi = .1;
Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
psihat = psi;

options = optimset('Display', 'off', 'LargeScale','off') ;

    % (legacy clock-seed line 38 deliberately dropped)
disp('Starting MCMC for BVAR-MA.... ');

for isim = 1:nsims + burnin

        %% sample Sig and A
    Xtld = Hpsi\X;
    Ytld = Hpsi\shortY;
    iO = sparse(1:T,1:T,[1/(1+psi^2) ones(1,T-1)]);
    XtldiO = Xtld'*iO;
    KA = sparse(1:k,1:k,1./VA0) + XtldiO*Xtld;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XtldiO*Ytld);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + Ytld'*iO*Ytld ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    Sig = iwishrnd(Shat,nu0+T);
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig';

    %% sample psi [legacy llike_MA -> bvt.ml.llike_ma]
    U = shortY - X*A;
    lp_psi = @(x) bvt.ml.llike_ma(x,U,Sig) + lpri_psi(x);
    if (mod(isim,100)==0) || isim == 1 %% get the Hessian every 100 iterations
        [psihat,fval,exitflag,output,grad,hess] ...
            = fminunc(@(x)-lp_psi(x),psihat,options);   %#ok<ASGLU>
        [tmpCpsi, flag] = chol(hess,'lower');           %#ok<ASGLU>
        if flag == 0
            Kpsic = hess;
        else
            Kpsic = 1/.05^2;
        end
    else
        psihat = fminbnd(@(x)-lp_psi(x),-.99,.99);
    end
    psic = psihat + 1/sqrt(Kpsic)*randn;
    if abs(psic)<.99
        alpMH =  lp_psi(psic) - lp_psi(psi) + ...
            -.5*(psi-psihat)^2*Kpsic + .5*(psic-psihat)^2*Kpsic;
    else
        alpMH = -inf;
    end
    if alpMH > log(rand)
        psi = psic;
        countpsi = countpsi + 1;
        Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
    end

    if isim > burnin
        isave = isim - burnin;
        store_A = store_A + A;
        store_Sig = store_Sig + Sig;
        store_psi(isave,:) = psi;
    end

    if ( mod(isim, pr.progress_every) ==0 )
        disp(  [ num2str(isim) ' loops... ' ] )
    end

end

    % script-tail posterior summaries [lines 112-114]
res = struct();
res.A_mean = store_A/nsims;
res.Sig_mean = store_Sig/nsims;
res.psi_mean = mean(store_psi)';
res.store_A = store_A; res.store_Sig = store_Sig;
res.store_psi = store_psi; res.store_ppsi = store_ppsi;
res.countpsi = countpsi;
res.pri = pri;
res.state = struct('psi', psi, 'psihat', psihat, 'Kpsic', Kpsic, ...
    'Sig', Sig, 'A', A, 'countpsi', countpsi);
end

% -------------------------------------------------------------------------
function res = mcmc_t_csv(shortY, X, pri, T, n, k, nsims, burnin, pr)
% BVAR_t_CSV.m functionized line-for-line (clock-seed line 43 dropped;
% figures not reproduced). Draw order per sweep:
% (Sig,A) -> h -> lam -> nu -> sigh2 -> rho.
A0 = pri.A0; VA0 = pri.VA0; nu0 = pri.nu0; S0 = pri.S0;
nuh0 = pr.nuh0; Sh0 = pr.Sh0;           % BVAR_t_CSV.m line 10
rho0 = pr.rho0; Vrho = pr.Vrho;         % line 11
nuub = pr.nuub;                         % line 12
pri.nuh0 = nuh0; pri.Sh0 = Sh0; pri.rho0 = rho0; pri.Vrho = Vrho; pri.nuub = nuub;

    % initialize for storage [lines 23-31]
store_Sig = zeros(n,n);
store_A = zeros(k,n);
store_h = zeros(nsims,T);
store_lam = zeros(nsims,T);
store_theta = zeros(nsims,3);

counth = 0; countrho = 0; countnu = 0;
nugrid = linspace(2,nuub,700)';
store_pnu = zeros(700,1);

    % initialize the chain [lines 34-39]
h = zeros(T,1);
nu = 5;
rho = .8;
sigh2 = .1;
Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T); %#ok<NASGU> % dead in the loop; kept verbatim
lam = 1./gamrnd(nu/2,2/nu,T,1);

    % (legacy clock-seed line 43 deliberately dropped)
disp('Starting MCMC for BVAR-t-CSV.... ');

for isim = 1:nsims + burnin

        %% sample Sig and A
    iOm = sparse(1:T,1:T,exp(-h)./lam);
    XiOm = X'*iOm;
    KA = sparse(1:k,1:k,1./VA0) + XiOm*X;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiOm*shortY);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + shortY'*iOm*shortY ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    Sig = iwishrnd(Shat,nu0+T);
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig';

        %% sample h [legacy sample_h -> bvt.sv.csv_armh]
    U = shortY - X*A;
    tmp = (U/CSig');
    s2 = sum(tmp.^2,2)./lam;
    [h, flag] = bvt.sv.csv_armh(s2,rho,sigh2,h,n);
    counth = counth + flag;

        %% sample lam
    U = shortY - X*A;
    tmp = (U/CSig');
    s2 = sum(tmp.^2,2)./exp(h);
    lam = 1./gamrnd((n+nu)/2,2./(s2+nu));

        %% sample nu [legacy sample_nu -> bvt.sv.nu_studentt]
    [nu, flag, fnu] = bvt.sv.nu_studentt(lam,nu,nuub);
    countnu = countnu + flag;

        %% sample sigh2
    eh = [h(1)*sqrt(1-rho^2);  h(2:end)-rho*h(1:end-1)];
    sigh2 = 1/gamrnd(nuh0+T/2,1/(Sh0 + sum(eh.^2)/2));

        %% sample rho [estimation truncation bound .9999 - see preset]
    Krho = 1/Vrho + sum(h(1:T-1).^2)/sigh2;
    rhohat = Krho\(rho0/Vrho + h(1:T-1)'*h(2:T)/sigh2);
    rhoc = rhohat + sqrt(Krho)'\randn;
    grho = @(x) -.5*log(sigh2./(1-x.^2))-.5*(1-x.^2)/sigh2*h(1)^2;
    if abs(rhoc)<.9999
        alpMH = exp(grho(rhoc)-grho(rho));
        if alpMH>rand
            rho = rhoc;
            countrho = countrho+1;
            Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T); %#ok<NASGU>
        end
    end


    if isim > burnin
        isave = isim - burnin;
        store_A = store_A + A;
        store_Sig = store_Sig + Sig;
        store_h(isave,:) = h';
        store_lam(isave,:) = lam';
        store_theta(isave,:) = [nu rho sigh2];

        % compute the density of nu
        tmpden = fnu(nugrid);
        tmpden = exp(tmpden-max(tmpden));
        tmpden = tmpden/(sum(tmpden)*(nugrid(2)-nugrid(1)));
        store_pnu = store_pnu + tmpden;
    end

    if ( mod(isim, pr.progress_every) ==0 )
        disp(  [ num2str(isim) ' loops... ' ] )
    end

end

    % script-tail posterior summaries [lines 121-126]
res = struct();
res.A_mean = store_A/nsims;
res.Sig_mean = store_Sig/nsims;
res.theta_mean = mean(store_theta)';
res.h_mean = mean(store_h)';
res.pnu_mean = store_pnu/nsims;
res.CSV_std_mean = mean(exp(store_h/2))';
res.store_A = store_A; res.store_Sig = store_Sig;
res.store_h = store_h; res.store_lam = store_lam;
res.store_theta = store_theta; res.store_pnu = store_pnu;
res.counth = counth; res.countrho = countrho; res.countnu = countnu;
res.pri = pri;
res.state = struct('h', h, 'lam', lam, 'nu', nu, 'rho', rho, 'sigh2', sigh2, ...
    'Sig', Sig, 'A', A, 'counth', counth, 'countrho', countrho, 'countnu', countnu);
end

% -------------------------------------------------------------------------
function res = mcmc_t_ma(shortY, X, pri, T, n, k, nsims, burnin, pr)
% BVAR_t_MA.m functionized line-for-line (clock-seed line 45 dropped;
% figures not reproduced). Draw order per sweep: (Sig,A) -> psi1 -> lam -> nu.
A0 = pri.A0; VA0 = pri.VA0; nu0 = pri.nu0; S0 = pri.S0;
psi0 = pr.psi0; Vpsi = pr.Vpsi;         % BVAR_t_MA.m line 9
lpri_psi = @(x) -.5*(x-psi0)^2/Vpsi -1e10*(x<-.99 || x>.99);    % line 10
nuub = pr.nuub;                         % line 12
pri.psi0 = psi0; pri.Vpsi = Vpsi; pri.nuub = nuub;

    % initialize for storage [lines 22-35]
store_Sig = zeros(n,n);
store_A = zeros(k,n);
store_lam = zeros(nsims,T);
store_theta = zeros(nsims,2);

options = optimset('Display', 'off', 'LargeScale','off') ;
countnu = 0; countpsi = 0;

ngrid = 300;
psigrid = linspace(-.1,.5,ngrid)';      %#ok<NASGU> % legacy line 32; read only by the commented-out ppsi block
store_ppsi = zeros(ngrid,1);
nugrid = linspace(2,50,ngrid)';         % NOTE: 2..50 here, not 2..nuub - legacy line 34
store_pnu = zeros(ngrid,1);

    % initialize the chain [lines 38-42]
nu = 5;
lam = 1./gamrnd(nu/2,2/nu,T,1);
psi1 = -.1;
Hpsi = speye(T) + psi1*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
psihat = psi1;

    % (legacy clock-seed line 45 deliberately dropped)
disp('Starting MCMC for BVAR-t-MA.... ');

for isim = 1:nsims + burnin

        %% sample Sig and A
    Xtld = Hpsi\X;
    Ytld = Hpsi\shortY;
    iO_lam = sparse(1:T,1:T,1./lam);
    XiO = Xtld'*iO_lam;
    KA = sparse(1:k,1:k,1./VA0) + XiO*Xtld;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiO*Ytld);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + Ytld'*iO_lam*Ytld ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    Sig = iwishrnd(Shat,nu0+T);
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig';

      %% sample psi1 [legacy llike_CSV_MA (ROOT copy) -> bvt.ml.llike_csv_ma,
      %  reused with h := log(lam)]
    U = shortY - X*A;
    lp_psi = @(x) bvt.ml.llike_csv_ma(x,U,Sig,log(lam)) + lpri_psi(x);
    if (mod(isim,100)==0) || isim == 1 %% get the Hessian every 100 iterations
        [psihat,fval,exitflag,output,grad,hess] ...
            = fminunc(@(x)-lp_psi(x),psihat,options);   %#ok<ASGLU>
        [tmpCpsi, flag] = chol(hess,'lower');           %#ok<ASGLU>
        if flag == 0
            Kpsic = hess;
        else
            Kpsic = 1/.05^2;
        end
    else
        psihat = fminbnd(@(x)-lp_psi(x),-.99,.99);
    end
    psic = psihat + 1/sqrt(Kpsic)*randn;
    if abs(psic)<.99
        alpMH =  lp_psi(psic) - lp_psi(psi1) + ...
            -.5*(psi1-psihat)^2*Kpsic + .5*(psic-psihat)^2*Kpsic;
    else
        alpMH = -inf;
    end
    if alpMH > log(rand)
        psi1 = psic;
        countpsi = countpsi + 1;
        Hpsi = speye(T) + psi1*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
    end

        %% sample lam
    Utld = Hpsi\(shortY - X*A);
    tmp = (Utld/CSig');
    s2 = sum(tmp.^2,2);
    s2(1) = s2(1)/(1+psi1^2);
    lam = 1./gamrnd((n+nu)/2,2./(s2+nu));

        %% sample nu [legacy sample_nu -> bvt.sv.nu_studentt]
    [nu, flag, fnu] = bvt.sv.nu_studentt(lam,nu,nuub);
    countnu = countnu + flag;

    if isim > burnin
        isave = isim - burnin;
        store_A = store_A + A;
        store_Sig = store_Sig + Sig;
        store_lam(isave,:) = lam';
        store_theta(isave,:) = [psi1 nu];

        % compute the posterior density of nu
        tmpden = fnu(nugrid);
        tmpden = exp(tmpden-max(tmpden));
        tmpden = tmpden/(sum(tmpden)*(nugrid(2)-nugrid(1)));
        store_pnu = store_pnu + tmpden;
    end

    if ( mod(isim, pr.progress_every) ==0 )
        disp(  [ num2str(isim) ' loops... ' ] )
    end

end

    % script-tail posterior summaries [lines 137-141]
res = struct();
res.A_mean = store_A/nsims;
res.Sig_mean = store_Sig/nsims;
res.theta_mean = mean(store_theta)';
res.ppsi_mean = store_ppsi/nsims;
res.pnu_mean = store_pnu/nsims;
res.store_A = store_A; res.store_Sig = store_Sig;
res.store_lam = store_lam; res.store_theta = store_theta;
res.store_ppsi = store_ppsi; res.store_pnu = store_pnu;
res.countnu = countnu; res.countpsi = countpsi;
res.pri = pri;
res.state = struct('psi1', psi1, 'psihat', psihat, 'Kpsic', Kpsic, ...
    'nu', nu, 'lam', lam, 'Sig', Sig, 'A', A, ...
    'countnu', countnu, 'countpsi', countpsi);
end

% -------------------------------------------------------------------------
function res = mcmc_csv_ma(shortY, X, pri, T, n, k, nsims, burnin, pr)
% BVAR_CSV_MA.m functionized line-for-line (clock-seed line 45 dropped;
% figures not reproduced). Draw order per sweep:
% (Sig,A) -> h -> sigh2 -> rho -> psi.
A0 = pri.A0; VA0 = pri.VA0; nu0 = pri.nu0; S0 = pri.S0;
psi0 = pr.psi0; Vpsi = pr.Vpsi;         % BVAR_CSV_MA.m line 9
lpri_psi = @(x) -.5*(x-psi0)^2/Vpsi -1e10*(x<-.99 || x>.99);    % line 10
nuh0 = pr.nuh0; Sh0 = pr.Sh0;           % line 12
rho0 = pr.rho0; Vrho = pr.Vrho;         % line 13
pri.psi0 = psi0; pri.Vpsi = Vpsi; pri.nuh0 = nuh0; pri.Sh0 = Sh0;
pri.rho0 = rho0; pri.Vrho = Vrho;

    % initialize for storage [lines 23-32]
store_Sig = zeros(n,n);
store_A = zeros(k,n);
store_h = zeros(nsims,T);
store_theta = zeros(nsims,3);
options = optimset('Display', 'off', 'LargeScale','off') ;
counth = 0; countrho = 0; countpsi = 0;
ngrid = 300;
psigrid = linspace(-.1,.4,ngrid)';      %#ok<NASGU> % legacy line 31; read only by the commented-out ppsi block
store_ppsi = zeros(ngrid,1);

    % initialize the chain [lines 34-42]
h = zeros(T,1);
psi = -.1;
rho = .8;
sigh2 = .1;
Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T); %#ok<NASGU> % dead in the loop; kept verbatim
Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);

psihat = psi;

    % (legacy clock-seed line 45 deliberately dropped)
disp('Starting MCMC for BVAR-CSV-MA.... ');

for isim = 1:nsims + burnin

        %% sample Sig and A
    Xtld = Hpsi\X;
    Ytld = Hpsi\shortY;
    iO_hpsi = sparse(1:T,1:T,[1/(1+psi^2)*exp(-h(1)); exp(-h(2:end))]);
    XiO = Xtld'*iO_hpsi;
    KA = sparse(1:k,1:k,1./VA0) + XiO*Xtld;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiO*Ytld);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + Ytld'*iO_hpsi*Ytld ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    Sig = iwishrnd(Shat,nu0+T);
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig';

        %% sample h [legacy sample_h -> bvt.sv.csv_armh]
    U = shortY - X*A;
    Utld = Hpsi\U;
    tmp = (Utld/CSig');
    s2 = sum(tmp.^2,2);
    s2(1) = s2(1)/(1+psi^2);
    [h, flag] = bvt.sv.csv_armh(s2,rho,sigh2,h,n);
    counth = counth + flag;

        %% sample sigh2
    eh = [h(1)*sqrt(1-rho^2);  h(2:end)-rho*h(1:end-1)];
    sigh2 = 1/gamrnd(nuh0+T/2,1/(Sh0 + sum(eh.^2)/2));

        %% sample rho [estimation truncation bound .9999 - see preset]
    Krho = 1/Vrho + sum(h(1:T-1).^2)/sigh2;
    rhohat = Krho\(rho0/Vrho + h(1:T-1)'*h(2:T)/sigh2);
    rhoc = rhohat + sqrt(Krho)'\randn;
    grho = @(x) -.5*log(sigh2./(1-x.^2))-.5*(1-x.^2)/sigh2*h(1)^2;
    if abs(rhoc)<.9999
        alpMH = exp(grho(rhoc)-grho(rho));
        if alpMH>rand
            rho = rhoc;
            countrho = countrho+1;
            Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T); %#ok<NASGU>
        end
    end

      %% sample psi [legacy llike_CSV_MA (ROOT copy) -> bvt.ml.llike_csv_ma]
    U = shortY - X*A;
    lp_psi = @(x) bvt.ml.llike_csv_ma(x,U,Sig,h) + lpri_psi(x);
    if (mod(isim,100)==0) || isim == 1 %% get the Hessian every 100 iterations
        [psihat,fval,exitflag,output,grad,hess] ...
            = fminunc(@(x)-lp_psi(x),psihat,options);   %#ok<ASGLU>
        [tmpCpsi, flag] = chol(hess,'lower');           %#ok<ASGLU>
        if flag == 0
            Kpsic = hess;
        else
            Kpsic = 1/.05^2;
        end
    else
        psihat = fminbnd(@(x)-lp_psi(x),-.99,.99);
    end
    psic = psihat + 1/sqrt(Kpsic)*randn;
    if abs(psic)<.99
        alpMH =  lp_psi(psic) - lp_psi(psi) + ...
            -.5*(psi-psihat)^2*Kpsic + .5*(psic-psihat)^2*Kpsic;
    else
        alpMH = -inf;
    end
    if alpMH > log(rand)
        psi = psic;
        countpsi = countpsi + 1;
        Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
    end

    if isim > burnin
        isave = isim - burnin;
        store_A = store_A + A;
        store_Sig = store_Sig + Sig;
        store_h(isave,:) = h';
        store_theta(isave,:) = [psi rho sigh2];
    end

    if ( mod(isim, pr.progress_every) ==0 )
        disp(  [ num2str(isim) ' loops... ' ] )
    end

end

    % script-tail posterior summaries [lines 147-152]
res = struct();
res.A_mean = store_A/nsims;
res.Sig_mean = store_Sig/nsims;
res.h_mean = mean(exp(store_h))';       % legacy line 149: mean of exp(h) here
res.theta_mean = mean(store_theta)';
res.ppsi_mean = store_ppsi/nsims;
res.CSV_std_mean = mean(exp(store_h/2))';
res.store_A = store_A; res.store_Sig = store_Sig;
res.store_h = store_h; res.store_theta = store_theta;
res.store_ppsi = store_ppsi;
res.counth = counth; res.countrho = countrho; res.countpsi = countpsi;
res.pri = pri;
res.state = struct('h', h, 'psi', psi, 'rho', rho, 'sigh2', sigh2, ...
    'psihat', psihat, 'Kpsic', Kpsic, 'Sig', Sig, 'A', A, ...
    'counth', counth, 'countrho', countrho, 'countpsi', countpsi);
end

% -------------------------------------------------------------------------
function res = mcmc_csv_t_ma(shortY, X, pri, T, n, k, nsims, burnin, pr)
% BVAR_CSV_t_MA.m functionized line-for-line (clock-seed line 51 dropped;
% figures not reproduced). Draw order per sweep:
% (Sig,A) -> lam -> h -> sigh2 -> rho -> psi -> nu.
% Verbatim legacy quirks kept (see preset and tests/variant_map.md): the rho
% MH truncation bound is .99 here (models 3/5/7 use .9999), and the lam step
% applies NO (1+psi^2) correction to the first observation's s2 (model 6
% does) - the h step does.
A0 = pri.A0; VA0 = pri.VA0; nu0 = pri.nu0; S0 = pri.S0;
psi0 = pr.psi0; Vpsi = pr.Vpsi;         % BVAR_CSV_t_MA.m line 9
lpri_psi = @(x) -.5*(x-psi0)^2/Vpsi -1e10*(x<-.99 || x>.99);    % line 10
nuh0 = pr.nuh0; Sh0 = pr.Sh0;           % line 12
rho0 = pr.rho0; Vrho = pr.Vrho;         % line 13
nuub = pr.nuub;                         % line 14
pri.psi0 = psi0; pri.Vpsi = Vpsi; pri.nuh0 = nuh0; pri.Sh0 = Sh0;
pri.rho0 = rho0; pri.Vrho = Vrho; pri.nuub = nuub;

    % initialize for storage [lines 25-36]
store_Sig = zeros(n,n);
store_A = zeros(k,n);
store_h = zeros(nsims,T);
store_lam = zeros(nsims,T);
store_theta = zeros(nsims,4);

options = optimset('Display', 'off', 'LargeScale','off') ;
ngrid = 300;
psigrid = linspace(-.1,.4,ngrid)';      %#ok<NASGU> % legacy line 33; read only by the commented-out ppsi block
store_ppsi = zeros(ngrid,1);
nugrid = linspace(2,nuub,ngrid)';
store_pnu = zeros(ngrid,1);

    % initialize the chain [lines 39-48]
h = zeros(T,1);
nu = 5;
psi = -.1;
rho = .8;
sigh2 = .1;
Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T); %#ok<NASGU> % dead in the loop; kept verbatim
Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
count_h = 0; count_rho = 0; count_psi = 0; count_nu = 0;
psihat = psi;
lam = 1./gamrnd(nu/2,2/nu,T,1);

    % (legacy clock-seed line 51 deliberately dropped)
disp('Starting MCMC for BVAR-CSV-t-MA.... ');

for isim = 1:nsims + burnin
        % sample Sig and A
    Xtld = Hpsi\X;
    Ytld = Hpsi\shortY;
    iO_h_lam_psi = sparse(1:T,1:T,[1/(1+psi^2)*exp(-h(1));exp(-h(2:end))]./lam);
    XiO = Xtld'*iO_h_lam_psi;
    KA = sparse(1:k,1:k,1./VA0) + XiO*Xtld;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiO*Ytld);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + Ytld'*iO_h_lam_psi*Ytld ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    Sig = iwishrnd(Shat,nu0+T);
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig';

        % sample lam [NO (1+psi^2) first-observation correction - verbatim]
    U = shortY - X*A;
    Utld = Hpsi\U;
    tmp = (Utld/CSig');
    s2_lam = sum(tmp.^2,2)./exp(h);
    lam = 1./gamrnd((n+nu)/2,2./(s2_lam+nu));

        % sample h [legacy sample_h -> bvt.sv.csv_armh; note: reuses `tmp`
        % from the lam step with the NEW lam - verbatim]
    s2_h = sum(tmp.^2,2)./lam;
    s2_h(1) = s2_h(1)/(1+psi^2);
    [h, flag] = bvt.sv.csv_armh(s2_h,rho,sigh2,h,n);
    count_h = count_h + flag;

        % sample sigh2
    eh = [h(1)*sqrt(1-rho^2);  h(2:end)-rho*h(1:end-1)];
    sigh2 = 1/gamrnd(nuh0+T/2,1/(Sh0 + sum(eh.^2)/2));

        % sample rho [estimation truncation bound .99 HERE - see preset]
    Krho = 1/Vrho + sum(h(1:T-1).^2)/sigh2;
    rhohat = Krho\(rho0/Vrho + h(1:T-1)'*h(2:T)/sigh2);
    rhoc = rhohat + sqrt(Krho)'\randn;
    grho = @(x) -.5*log(sigh2./(1-x.^2))-.5*(1-x.^2)/sigh2*h(1)^2;
    if abs(rhoc)<.99
        alpMH = exp(grho(rhoc)-grho(rho));
        if alpMH>rand
            rho = rhoc;
            count_rho = count_rho+1;
            Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T); %#ok<NASGU>
        end
    end

        % sample psi [legacy llike_CSV_MA (ROOT copy) -> bvt.ml.llike_csv_ma]
    U_psi = U./repmat(sqrt(lam),1,n);
    lp_psi = @(x) bvt.ml.llike_csv_ma(x,U_psi,Sig,h) + lpri_psi(x);
    if (mod(isim,100)==0) || isim == 1 %% get the Hessian every 100 iterations
        [psihat,fval,exitflag,output,grad,hess] ...
            = fminunc(@(x)-lp_psi(x),psihat,options);   %#ok<ASGLU>
        [tmpCpsi, flag] = chol(hess,'lower');           %#ok<ASGLU>
        if flag == 0
            Kpsic = hess;
        else
            Kpsic = 1/.05^2;
        end
    else
        psihat = fminbnd(@(x)-lp_psi(x),-.99,.99);
    end
    psic = psihat + 1/sqrt(Kpsic)*randn;
    if abs(psic)<.99
        alpMH =  lp_psi(psic) - lp_psi(psi) + ...
            -.5*(psi-psihat)^2*Kpsic + .5*(psic-psihat)^2*Kpsic;
    else
        alpMH = -inf;
    end
    if alpMH > log(rand)
        psi = psic;
        count_psi = count_psi + 1;
        Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
    end

        % sample nu [legacy sample_nu -> bvt.sv.nu_studentt]
    [nu,flag,fnu] = bvt.sv.nu_studentt(lam,nu,nuub);
    count_nu = count_nu + flag;

    if isim > burnin
        isave = isim - burnin;
        store_A = store_A + A;
        store_Sig = store_Sig + Sig;
        store_h(isave,:) = h';
        store_lam(isave,:) = lam';
        store_theta(isave,:) = [psi rho sigh2 nu];

        % compute the density of nu
        tmpden = fnu(nugrid);
        tmpden = exp(tmpden-max(tmpden));
        tmpden = tmpden/(sum(tmpden)*(nugrid(2)-nugrid(1)));
        store_pnu = store_pnu + tmpden;
    end

    if ( mod(isim, pr.progress_every) ==0 )
        disp(  [ num2str(isim) ' loops... ' ] )
    end

end

    % script-tail posterior summaries [lines 169-174]
res = struct();
res.A_mean = store_A/nsims;
res.Sig_mean = store_Sig/nsims;
res.CSV_std_mean = mean(exp(store_h/2))';
res.theta_mean = mean(store_theta)';
res.ppsi_mean = store_ppsi/nsims;
res.pnu_mean = store_pnu/nsims;
res.store_A = store_A; res.store_Sig = store_Sig;
res.store_h = store_h; res.store_lam = store_lam;
res.store_theta = store_theta;
res.store_ppsi = store_ppsi; res.store_pnu = store_pnu;
res.count_h = count_h; res.count_rho = count_rho;
res.count_psi = count_psi; res.count_nu = count_nu;
res.pri = pri;
res.state = struct('h', h, 'lam', lam, 'psi', psi, 'rho', rho, ...
    'sigh2', sigh2, 'nu', nu, 'psihat', psihat, 'Kpsic', Kpsic, ...
    'Sig', Sig, 'A', A, 'count_h', count_h, 'count_rho', count_rho, ...
    'count_psi', count_psi, 'count_nu', count_nu);
end
