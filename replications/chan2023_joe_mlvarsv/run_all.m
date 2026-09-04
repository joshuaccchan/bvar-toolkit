% chan2023_joe_mlvarsv/run_all - functionized estimation pipeline of Chan
% (2023, JoE): the five VARs whose marginal likelihoods the paper compares -
% homoskedastic natural-conjugate (VAR-NCP), common SV (VAR-CSV), Cholesky SV
% (VAR-SV), factor SV (VAR-FSV), and Cholesky SV with an outlier component
% (VAR-SVO).
%
%   out = run_all(model, is_kappafixed, is_kappasym, nsim, burnin, seed, varid)
%
%   model         - 'VAR-NCP' | 'VAR-CSV' | 'VAR-SV' | 'VAR-FSV' | 'VAR-SVO',
%                   or the legacy numeric code 1-5 (main_varsv.m line 20);
%                   default 'VAR-CSV' (the legacy default model = 2)
%   is_kappafixed - hold the shrinkage hyperparameters fixed (main_varsv.m 22).
%                   Note VAR-SV/VAR-SVO still draw kappa4 - the legacy block
%                   sits before the branch (preset pr.sv.kappa4_always_drawn)
%   is_kappasym   - kappa2 = kappa1 (main_varsv.m 23); ignored by VAR-NCP and,
%                   as in the legacy switch, by VAR-CSV
%   nsim, burnin  - defaults 10000 / 1000 (main_varsv.m 26-27)
%   seed          - when nonempty, rng(seed,'twister') is set before any draw;
%                   omitted/empty uses the ambient stream as-is
%   varid         - data columns; default the active n = 15 selection
%                   (main_varsv.m 35). preset also carries the commented n = 7
%                   and n = 30 selections.
%
% Functionized 2026-09-03 (step 9). Reproduces main_varsv.m -> VAR_NCP.m /
% VAR_CSV.m / VAR_ARSV_redu.m / VAR_FSV.m / VAR_ARSVO_redu.m draw-for-draw
% bitwise (tests/unit/test_mlvarsv_equivalence.m) with one deliberate
% divergence: the four MCMC scripts re-seed the global stream from the wall
% clock (VAR_CSV.m 32, VAR_ARSV_redu.m 38, VAR_FSV.m 34, VAR_ARSVO_redu.m 45),
% which is irreproducible and switches MATLAB to the legacy v4/v5 generators;
% run_all drops those lines so the caller controls seeding. Every rng draw sits
% after that point, so one seed aligns a whole run, chain init included. The
% wall-clock timing displays and VAR_CSV.m's exp(h/2) figure are not reproduced.
%
% Scope: estimation only. The marginal-likelihood routines utility/ml_var_*.m
% are a separate phase and are not called here. Model 1 is the exception - its
% log-ML is inline and analytic in VAR_NCP.m, so run_all returns it.
%
% Core used: bvar.priors.niw('mlvarsv_ncp') (legacy prior_NCP), bvar.priors.minn
% (prior_Minn, n0pre = 4), bvar.priors.impact_B0 (prior_B0), bvar.priors.minnesota_C
% (get_C), bvar.sv.init_approx1N (getARh_approx1N), bvar.sv.csv_armh (sample_CSV),
% bvar.sv.ksc_ar1_mean (sample_SV), bvar.sv.sv_params (sample_SVpara, called at
% this package's phi bound .998), bvar.sv.svo_outlier (the o/po block),
% bvar.samplers.eq_var_redu_tri (the reduced-form coefficient draw, with the
% outlier scaling for VAR-SVO), bvar.samplers.alp_tri_cs (the B0 block; returns a
% row, transposed here to the legacy column), bvar.samplers.factor_fsv and
% bvar.samplers.eq_fsv_load (the FSV factor and loading blocks),
% bvar.util.build_lags / vec / vech / ldet / mgammaln, third_party/gigrnd.
% The VAR-CSV (Sig,A) natural-conjugate joint draw has no core counterpart yet
% - it is the same block chan2020_jbes_kronecker/run_all.m keeps inline, held
% for the springer family pass - so it stays verbatim in mcmc_csv below.
%
% See:
% Chan, J.C.C. (2023). Comparing stochastic volatility specifications for large
% Bayesian VARs, Journal of Econometrics, 235(2), 1419-1446.

function out = run_all(model, is_kappafixed, is_kappasym, nsim, burnin, seed, varid)
thisdir = fileparts(mfilename('fullpath'));

    % make bvar.* and gigrnd resolvable when called standalone
if isempty(which('bvar.samplers.eq_var_redu_tri'))
    root = fileparts(fileparts(thisdir));
    addpath(fullfile(root, 'core'));
end
if isempty(which('gigrnd'))
    root = fileparts(fileparts(thisdir));
    addpath(fullfile(root, 'third_party'));
end

    % constants: preset.m in this folder (cd guard pins name resolution)
od = cd(thisdir);
guard = onCleanup(@() cd(od));
pr = preset();
clear guard

names = {'VAR-NCP', 'VAR-CSV', 'VAR-SV', 'VAR-FSV', 'VAR-SVO'};
if nargin < 1 || isempty(model), model = 'VAR-CSV'; end
if isnumeric(model)
    assert(any(model == 1:5), 'run_all:model', 'numeric model must be 1-5');
    imodel = model;
else
    imodel = find(strcmpi(model, names), 1);
    if isempty(imodel)
        error('run_all:model', 'model must be one of %s (or 1-5)', strjoin(names, ', '));
    end
end
if nargin < 2 || isempty(is_kappafixed), is_kappafixed = false;          end
if nargin < 3 || isempty(is_kappasym),   is_kappasym   = false;          end
if nargin < 4 || isempty(nsim),          nsim   = pr.nsim_default;       end
if nargin < 5 || isempty(burnin),        burnin = pr.burnin_default;     end
if nargin < 6, seed = []; end
if nargin < 7 || isempty(varid),         varid = pr.varid;               end
if ~isempty(seed)
    rng(seed, 'twister');
end

    % data and design [main_varsv.m 33-53]
p = pr.p;
r = pr.r;
data_all = load(fullfile(thisdir, 'legacy', pr.data_file));   % legacy folder, read-only
data = data_all(:, varid);
Y0 = data(1:pr.n0, :);
Y  = data(pr.n0+1:end, :);
[T, n] = size(Y);
[~, X] = bvar.util.build_lags([Y0(end-p+1:end, :); Y], p);   % identical to the legacy inline loop
k_beta = n*(n-1)/2;         % dimension of B0
k_alp  = n^2*p + n;         % dimension of A
k = k_alp/n;
    % legacy line 44 also builds y = reshape(Y',n*T,1); no estimation script
    % reads it (grep-verified) - not reproduced.

    % priors and model name [main_varsv.m 56-122]
kappa = pr.kappa; kappa1 = pr.kappa1; kappa3 = pr.kappa3; kappa4 = pr.kappa4;
if is_kappasym
    kappa2 = kappa1;
else
    kappa2 = pr.kappa2;
end
sig2_hat = [];
Hyper = struct();
switch imodel
    case 1
        model_name = 'VAR-NCP';
        [Hyper.A0, Hyper.VA, Hyper.nu0, Hyper.S0] = ...
            bvar.priors.niw(p, [kappa kappa3], Y0, Y, 'mlvarsv_ncp');
    case 2
        if is_kappafixed
            model_name = 'VAR-CSV-fixed-kappa';
        else
            model_name = 'VAR-CSV';
        end
        [Hyper.A0, Hyper.VA, Hyper.nu0, Hyper.S0] = ...
            bvar.priors.niw(p, [kappa kappa3], Y0, Y, 'mlvarsv_ncp');
        Hyper.c0 = pr.csv.c0;
        Hyper.nuh = pr.nuh_scalar; Hyper.Sh = pr.Sh_factor*(Hyper.nuh-1);
        Hyper.phi0 = pr.phi0_scalar; Hyper.Vphi = pr.Vphi_scalar;
    case 3
        if is_kappafixed
            model_name = 'VAR-SV-fixed-kappa';
        elseif is_kappasym
            model_name = 'VAR-SV-sym';
        else
            model_name = 'VAR-SV';
        end
        [Hyper.alp0, Hyper.Valp, sig2_hat] = ...
            bvar.priors.minn(p, kappa1, kappa2, kappa3, Y0, Y, 4);
        [Hyper.beta0, Hyper.Vbeta] = bvar.priors.impact_B0(Y0, Y, kappa4);
        Hyper.c0 = pr.sv.c0;
        Hyper.nuh = pr.nuh_scalar*ones(n,1); Hyper.Sh = pr.Sh_factor*(Hyper.nuh-1);
        Hyper.mu0 = pr.mu0_scalar*ones(n,1); Hyper.Vmu = pr.Vmu_scalar*ones(n,1);
        Hyper.phi0 = pr.phi0_scalar*ones(n,1); Hyper.Vphi = pr.Vphi_scalar*ones(n,1);
    case 4
        if is_kappafixed
            model_name = 'VAR-FSV-fixed-kappa';
        elseif is_kappasym
            model_name = 'VAR-FSV-sym';
        else
            model_name = 'VAR-FSV';
        end
        [Hyper.alp0, Hyper.Valp, sig2_hat] = ...
            bvar.priors.minn(p, kappa1, kappa2, kappa3, Y0, Y, 4);
        Hyper.c0 = pr.fsv.c0;
        Hyper.nuh = pr.nuh_scalar*ones(n+r,1); Hyper.Sh = pr.Sh_factor*(Hyper.nuh-1);
        Hyper.mu0 = pr.mu0_scalar*ones(n+r,1); Hyper.Vmu = pr.Vmu_scalar*ones(n+r,1);
        Hyper.phi0 = pr.phi0_scalar*ones(n+r,1); Hyper.Vphi = pr.Vphi_scalar*ones(n+r,1);
        Hyper.l0 = pr.fsv.l0;
        Hyper.Vl = pr.fsv.Vl;
    case 5
        if is_kappafixed
            model_name = 'VAR-SVO-fixed-kappa';
        elseif is_kappasym
            model_name = 'VAR-SVO-sym';
        else
            model_name = 'VAR-SVO';
        end
        [Hyper.alp0, Hyper.Valp, sig2_hat] = ...
            bvar.priors.minn(p, kappa1, kappa2, kappa3, Y0, Y, 4);
        [Hyper.beta0, Hyper.Vbeta] = bvar.priors.impact_B0(Y0, Y, kappa4);
        Hyper.c0 = pr.svo.c0;
        Hyper.nuh = pr.nuh_scalar*ones(n,1); Hyper.Sh = pr.Sh_factor*(Hyper.nuh-1);
        Hyper.mu0 = pr.mu0_scalar*ones(n,1); Hyper.Vmu = pr.Vmu_scalar*ones(n,1);
        Hyper.phi0 = pr.phi0_scalar*ones(n,1); Hyper.Vphi = pr.Vphi_scalar*ones(n,1);
        Hyper.p0a = pr.svo.p0a; Hyper.p0b = pr.svo.p0b;
end

disp(['Starting MCMC for ' model_name '...']);   % main_varsv.m 125
dims = struct('T',T, 'n',n, 'p',p, 'r',r, 'k',k, 'k_alp',k_alp, 'k_beta',k_beta);
switch imodel
    case 1
        res = ncp(Y, X, dims, Hyper);
    case 2
        res = mcmc_csv(Y, X, Y0, dims, Hyper, kappa, kappa3, is_kappafixed, nsim, burnin, pr);
    case 3
        res = mcmc_arsv(Y, X, Y0, dims, Hyper, sig2_hat, kappa1, kappa2, kappa3, kappa4, ...
            is_kappafixed, is_kappasym, nsim, burnin, pr);
    case 4
        res = mcmc_fsv(Y, X, Y0, dims, Hyper, sig2_hat, kappa1, kappa2, kappa3, ...
            is_kappafixed, is_kappasym, nsim, burnin, pr);
    case 5
        res = mcmc_arsvo(Y, X, Y0, dims, Hyper, sig2_hat, kappa1, kappa2, kappa3, kappa4, ...
            is_kappafixed, is_kappasym, nsim, burnin, pr);
end

out = res;
out.model = names{imodel};
out.model_num = imodel;
out.model_name = model_name;
out.is_kappafixed = is_kappafixed;
out.is_kappasym = is_kappasym;
out.nsim = nsim; out.burnin = burnin; out.seed = seed;
out.varid = varid;
out.T = T; out.n = n; out.p = p; out.r = r;
out.k = k; out.k_alp = k_alp; out.k_beta = k_beta;
out.Y = Y; out.X = X; out.Y0 = Y0;   % the ML phase needs them (run_ml.m)
out.Hyper = Hyper;
out.preset = pr;
end

% -------------------------------------------------------------------------
function res = ncp(Y, X, dims, Hyper)
% VAR_NCP.m functionized line-for-line. No MCMC, no rng. The cp_ml block
% (lines 18-21) is inline and analytic, so it is reproduced here.
T = dims.T; n = dims.n; k = dims.k;

XX = X'*X;                                              % VAR_NCP.m 10
A_tilde = XX\(X'*Y);                                    % 11
K_A = sparse(1:k,1:k,1./Hyper.VA) + XX;                 % 12
    % posterior mean of the VAR coefficients, arranged as a k by n matrix
A_hat = K_A\(sparse(1:k,1:k,Hyper.VA)\Hyper.A0 + XX*A_tilde);   % 14
S_hat = Hyper.S0 + Hyper.A0'*sparse(1:k,1:k,1./Hyper.VA)*Hyper.A0 + Y'*Y - A_hat'*K_A*A_hat;  % 15
S_hat = (S_hat+S_hat')/2;                               % 16

lml = -n*T/2*log(pi) - n/2*(sum(log(Hyper.VA)) + bvar.util.ldet(K_A)) ...
    + Hyper.nu0/2*bvar.util.ldet(Hyper.S0) ...
    - (Hyper.nu0+T)/2*bvar.util.ldet(S_hat) ...
    + bvar.util.mgammaln(n,(Hyper.nu0+T)/2) - bvar.util.mgammaln(n,Hyper.nu0/2);  % 19-20

res = struct('A_tilde',A_tilde, 'K_A',K_A, 'A_hat',A_hat, 'S_hat',S_hat, 'lml',lml);
end

% -------------------------------------------------------------------------
function res = mcmc_csv(Y, X, Y0, dims, Hyper, kappa, kappa3, is_kappafixed, nsim, burnin, pr)
% VAR_CSV.m functionized line-for-line (clock-seed line 32 dropped). Draw order
% per sweep: (Sig,A) -> h -> (phi,sig2) -> kappa.
T = dims.T; n = dims.n; p = dims.p; k = dims.k;

    % storage [VAR_CSV.m 9-13]
store_Sig = zeros(nsim,n*(n+1)/2);
store_a = zeros(nsim,k*n);
store_h = zeros(nsim,T);
store_hpara = zeros(nsim,2);
store_kappa = zeros(nsim,1);

    % initialize the chain [16-29]
sig2 = 1./gamrnd(Hyper.nuh,1./Hyper.Sh);
phi = Hyper.phi0;
[Hyper.A0,Hyper.VA] = bvar.priors.niw(p,[kappa kappa3],Y0,Y,'mlvarsv_ncp');
C = Hyper.VA/kappa; % the constant part of VA
K_A = sparse(1:k,1:k,1./Hyper.VA) + X'*X;
A = K_A\(sparse(1:k,1:k,Hyper.VA)\Hyper.A0 + X'*Y);
Sig = (Hyper.S0 + Hyper.A0'*sparse(1:k,1:k,1./Hyper.VA)*Hyper.A0 ...
        + Y'*Y - A'*K_A*A)/T;
U = Y - X*A;
tmp = (U/chol(Sig,'lower')');
s2 = sum(tmp.^2,2);
h = bvar.sv.csv_armh(s2,phi,sig2,zeros(T,1),n,true);
Hphi = speye(T) - phi*sparse(2:T,1:(T-1),ones(1,T-1),T,T); %#ok<NASGU> % dead in the legacy too [28]
count_h = 0; count_phi = 0;

for isim = 1:nsim + burnin
        % sample Sig and A [36-46]
    [~,Hyper.VA] = bvar.priors.niw(p,[kappa kappa3],Y0,Y,'mlvarsv_ncp');
    iOh = sparse(1:T,1:T,exp(-h));
    XiOh = X'*iOh;
    K_A = sparse(1:k,1:k,1./Hyper.VA) + XiOh*X;
    A_hat = K_A\(sparse(1:k,1:k,Hyper.VA)\Hyper.A0 + XiOh*Y);
    S_hat = Hyper.S0 + Hyper.A0'*sparse(1:k,1:k,1./Hyper.VA)*Hyper.A0 ...
        + Y'*iOh*Y - A_hat'*K_A*A_hat;
    S_hat = (S_hat+S_hat')/2; % adjust for rounding errors
    Sig = iwishrnd(S_hat,Hyper.nu0+T);
    CSig = chol(Sig,'lower');
    A = A_hat + (chol(K_A,'lower')'\randn(k,n))*CSig';

        % sample h [49-57 -> bvar.sv.csv_armh]
    U = Y - X*A;
    tmp = (U/CSig');
    s2 = sum(tmp.^2,2);
    if isim < pr.csv.forced_accept_before
        h = bvar.sv.csv_armh(s2,phi,sig2,h,n,true);
        flag = 0;
    else
        [h,flag] = bvar.sv.csv_armh(s2,phi,sig2,h,n);
    end
    count_h = count_h + flag;

        % sample phi and sig2 [61-63; mu passed as 0 gates the mu block off]
    [~,phi,sig2,flag_phi] = bvar.sv.sv_params(h,0,phi,Hyper,pr.phi_mh_bnd);
    count_phi = count_phi + flag_phi;
    Hphi = speye(T) - phi*sparse(2:T,1:(T-1),ones(1,T-1),T,T); %#ok<NASGU> % dead in the legacy too [63]

        % sample kappa [66-70]
    if ~is_kappafixed
        Q = diag((A-Hyper.A0)*(Sig\(A-Hyper.A0)'));
        tmpc = sum(Q(2:end)./C(2:end));
        kappa = gigrnd(Hyper.c0(1)-n^2*p/2,2*Hyper.c0(2),tmpc,1); % k=np+1
    end

    if isim > burnin                                    % [72-79]
        isave = isim - burnin;
        store_a(isave,:,:) = A(:);
        store_Sig(isave,:) = bvar.util.vech(Sig);
        store_h(isave,:) = h';
        store_hpara(isave,:) = [phi sig2];
        store_kappa(isave,:) = kappa;
    end

    if ( mod(isim, pr.progress_every) ==0 )             % [81-83]
        disp(  [ num2str(isim) ' loops... ' ] )
    end
end

    % posterior summaries [96-100]
res = struct();
res.store_Sig = store_Sig;
res.store_a = store_a;
res.store_h = store_h;
res.store_hpara = store_hpara;
res.store_kappa = store_kappa;
res.count_h = count_h;
res.count_phi = count_phi;
res.h_mean = mean(store_h)';
res.hpara_mean = mean(store_hpara)';
res.kappa_mean = mean(store_kappa);
res.CSV_std_mean = mean(exp(store_h/2))';
end

% -------------------------------------------------------------------------
function res = mcmc_arsv(Y, X, Y0, dims, Hyper, sig2_hat, kappa1, kappa2, kappa3, kappa4, ...
    is_kappafixed, is_kappasym, nsim, burnin, pr)
% VAR_ARSV_redu.m functionized line-for-line (clock-seed line 38 dropped). Draw
% order per sweep: alp -> beta -> h -> (mu,phi,sig2) -> kappa4, kappa1, kappa2.
T = dims.T; n = dims.n; p = dims.p; k = dims.k; k_alp = dims.k_alp; k_beta = dims.k_beta;

    % storage [VAR_ARSV_redu.m 9-14]
store_alp = zeros(nsim,k_alp);
store_beta = zeros(nsim,k_beta);
store_h = zeros(nsim,T,n);
store_hpara = zeros(nsim,3*n);
store_kappa = zeros(nsim,3);
count_phi = zeros(n,1);

    % initialize the Markov chain [17-35]
mu = log(sig2_hat);
sig2 = 1./gamrnd(Hyper.nuh,1./Hyper.Sh);
phi = Hyper.phi0;
A = zeros(k,n);
h = zeros(T,n);
XX = X'*X;
for ii=1:n
    iValpi = sparse(1:k,1:k,1./Hyper.Valp((ii-1)*k+1:ii*k));
    A(:,ii) = (XX + iValpi)\(X'*Y(:,ii));
    s2i = (Y(:,ii) - X*A(:,ii)).^2;
    mu(ii) = mean(log(s2i));
    h(:,ii) = bvar.sv.init_approx1N(s2i,mu(ii),phi(ii),sig2(ii));
end
alp = reshape(A,k_alp,1);                               %#ok<NASGU> % recomputed every sweep [30]
beta = zeros(k_beta,1);                                 %#ok<NASGU> % overwritten every sweep [31]
B0_id = nonzeros(tril(reshape(1:n^2,n,n),-1)');
B0 = eye(n);
[C_alp,idx_kappa1,idx_kappa2] = bvar.priors.minnesota_C(n,p,sig2_hat);
C_beta = Hyper.Vbeta/kappa4; % the constant part of Vbeta

for isim = 1:nsim + burnin
        % sample alp [42-58 -> bvar.samplers.eq_var_redu_tri]
    [~,Hyper.Valp] = bvar.priors.minn(p,kappa1,kappa2,kappa3,Y0,Y,4);
    U = zeros(T,n);                                     %#ok<NASGU> % dead in the legacy too [43]
    A = bvar.samplers.eq_var_redu_tri(Y,X,B0,h,A,Hyper.Valp,Hyper.alp0);
    alp = reshape(A,k_alp,1);

        % sample beta [61-74 -> bvar.samplers.alp_tri_cs; transposed to the
        % legacy column orientation, which store_beta and tmpc4 depend on]
    [~,Hyper.Vbeta] = bvar.priors.impact_B0(Y0,Y,kappa4);
    E = Y - X*A;
    beta = bvar.samplers.alp_tri_cs(E,h,Hyper.Vbeta)';
    B0(B0_id) = beta;

        % sample h [77-81]
    B0E = E*B0';
    for ii=1:n
        ystar = log(B0E(:,ii).^2 + pr.sv_offset);
        h(:,ii) = bvar.sv.ksc_ar1_mean(ystar,h(:,ii),mu(ii),phi(ii),sig2(ii));
    end

        % sample mu, phi and sig2 [84-85; ml_varsv phi bound .998]
    [mu,phi,sig2,flag_phi] = bvar.sv.sv_params(h,mu,phi,Hyper,pr.phi_mh_bnd);
    count_phi = count_phi + flag_phi;

        % sample kappa1, kappa2 and kappa4 [88-102]
    tmpc4 = sum(beta.^2./C_beta);
    kappa4 = gigrnd(Hyper.c0(3,1)-n*(n-1)/4,2*Hyper.c0(3,2),tmpc4,1);
    if is_kappafixed
        % nothing to do here
    elseif is_kappasym
        tmpc1 = sum(alp(idx_kappa1).^2./C_alp(idx_kappa1)) + ...
            sum(alp(idx_kappa2).^2./C_alp(idx_kappa2));
        kappa1 = gigrnd(Hyper.c0(1,1)-n^2*p/2,2*Hyper.c0(1,2),tmpc1,1);
        kappa2 = kappa1;
    else
        tmpc1 = sum(alp(idx_kappa1).^2./C_alp(idx_kappa1));
        tmpc2 = sum(alp(idx_kappa2).^2./C_alp(idx_kappa2));
        kappa1 = gigrnd(Hyper.c0(1,1)-n*p/2,2*Hyper.c0(1,2),tmpc1,1);
        kappa2 = gigrnd(Hyper.c0(2,1)-(n-1)*n*p/2,2*Hyper.c0(2,2),tmpc2,1);
    end

    if ( mod(isim, pr.progress_every) ==0 )             % [104-106]
        disp(  [ num2str(isim) ' loops... ' ] )
    end

    if isim > burnin                                    % [108-115]
        isave = isim-burnin;
        store_alp(isave,:) = reshape(A,1,k_alp);
        store_beta(isave,:) = beta';
        store_h(isave,:,:) = h;
        store_hpara(isave,:) = [mu',phi',sig2'];
        store_kappa(isave,:) = [kappa1,kappa2,kappa4];
    end
end

    % posterior summaries [119-124]
res = struct();
res.store_alp = store_alp;
res.store_beta = store_beta;
res.store_h = store_h;
res.store_hpara = store_hpara;
res.store_kappa = store_kappa;
res.count_phi = count_phi;
res.alp_mean = mean(store_alp)';
res.A_mean = reshape(res.alp_mean,k,n)';
res.beta_mean = mean(store_beta)';
res.h_mean = squeeze(mean(store_h));
res.hpara_mean = mean(store_hpara)';
res.kappa_mean = mean(store_kappa)';
end

% -------------------------------------------------------------------------
function res = mcmc_fsv(Y, X, Y0, dims, Hyper, sig2_hat, kappa1, kappa2, kappa3, ...
    is_kappafixed, is_kappasym, nsim, burnin, pr)
% VAR_FSV.m functionized line-for-line (clock-seed line 34 dropped). Draw order
% per sweep: F -> (A,L) -> h -> (mu,phi,sig2) -> kappa1, kappa2.
T = dims.T; n = dims.n; p = dims.p; r = dims.r; k = dims.k;

    % storage [VAR_FSV.m 9-16]
store_l = zeros(nsim,n*r-r*(r+1)/2);
store_A = zeros(k,n);
store_h = zeros(nsim,T,n+r);
store_F = zeros(nsim,T,r);
store_hpara = zeros(nsim,3*(n+r)); % [mu' phi sig2];
count_phi = zeros(n+r,1);
store_kappa = zeros(nsim,2);
L_idx = find(tril(ones(n,r),-1)~=0); % index of free elements of L

    % initialize the Markov chain [19-31]
A = (X'*X + sparse(1:k,1:k,1./[pr.fsv.A_init_VA(1),pr.fsv.A_init_VA(2)*ones(1,k-1)]))\(X'*Y);
L = [eye(r); pr.fsv.L_init_lower*ones(n-r,r)];
E = Y-X*A;                                              %#ok<NASGU> % dead in the legacy too [21]
varY = var(Y)';
mu = [log(varY/2); log(mean(varY))*ones(r,1)];
phi = pr.fsv.phi_init(1) + pr.fsv.phi_init(2)*rand(n+r,1);
sig2 = pr.fsv.sig2_init(1) + pr.fsv.sig2_init(2)*rand(n+r,1);
Uy = Y-X*A;
h = [zeros(T,n) repmat(mu(n+1:end)',T,1)];
for ii=1:n
    h(:,ii) = bvar.sv.init_approx1N(Uy(:,ii).^2,mu(ii),phi(ii),sig2(ii));
end
[C_alp,idx_kappa1,idx_kappa2] = bvar.priors.minnesota_C(n,p,sig2_hat);

for isim = 1:nsim + burnin
        % sample F [38-44 -> bvar.samplers.factor_fsv]
    F = bvar.samplers.factor_fsv(Y,X,A,L,h);

        % sample L and A [47-71 -> bvar.samplers.eq_fsv_load]
    [~,Hyper.Valp] = bvar.priors.minn(p,kappa1,kappa2,kappa3,Y0,Y,4);
    [A,L] = bvar.samplers.eq_fsv_load(Y,X,F,h,A,L,Hyper.Valp,Hyper.alp0,Hyper.Vl,Hyper.l0);
    alp = A(:);

        % sample h [75-79]
    Uy = Y -X*A -F*L';
    ystar = log([Uy F].^2 + pr.sv_offset);
    for jj=1:n+r
        h(:,jj) = bvar.sv.ksc_ar1_mean(ystar(:,jj),h(:,jj),mu(jj),phi(jj),sig2(jj));
    end

        % sample mu, phi and sig2 [82-83; ml_varsv phi bound .998]
    [mu,phi,sig2,flag_phi] = bvar.sv.sv_params(h,mu,phi,Hyper,pr.phi_mh_bnd);
    count_phi = count_phi + flag_phi;

        % sample kappa1 and kappa2 [86-98]
    if is_kappafixed
        % do nothing
    elseif is_kappasym
        tmpc1 = sum(alp(idx_kappa1).^2./C_alp(idx_kappa1)) + ...
            sum(alp(idx_kappa2).^2./C_alp(idx_kappa2));
        kappa1 = gigrnd(Hyper.c0(1,1)-n^2*p/2,2*Hyper.c0(1,2),tmpc1,1);
        kappa2 = kappa1;
    else
        tmpc1 = sum(alp(idx_kappa1).^2./C_alp(idx_kappa1));
        tmpc2 = sum(alp(idx_kappa2).^2./C_alp(idx_kappa2));
        kappa1 = gigrnd(Hyper.c0(1,1)-n*p/2,2*Hyper.c0(1,2),tmpc1,1);
        kappa2 = gigrnd(Hyper.c0(2,1)-(n-1)*n*p/2,2*Hyper.c0(2,2),tmpc2,1);
    end

    if isim > burnin                                    % [100-108]
        isave = isim - burnin;
        store_l(isave,:) = L(L_idx);
        store_A = store_A + A;
        store_h(isave,:,:) = h;
        store_F(isave,:,:) = F;
        store_hpara(isave,:) =  [mu',phi',sig2'];
        store_kappa(isave,:) = [kappa1,kappa2];
    end

    if ( mod(isim, pr.progress_every) ==0 )             % [110-112]
        disp(  [ num2str(isim) ' loops... ' ] )
    end
end

    % posterior summaries [116-122]
res = struct();
res.store_l = store_l;
res.store_A = store_A;
res.store_h = store_h;
res.store_F = store_F;
res.store_hpara = store_hpara;
res.store_kappa = store_kappa;
res.count_phi = count_phi;
res.l_mean = mean(store_l);
L_mean = [eye(r);ones(n-r,r)]; L_mean(L_idx) = res.l_mean;
res.L_mean = L_mean;
res.A_mean = store_A/nsim;
res.h_mean = squeeze(mean(store_h));
res.F_mean = squeeze(mean(store_F));
res.hpara_mean = mean(store_hpara)';
res.kappa_mean = mean(store_kappa)';
end

% -------------------------------------------------------------------------
function res = mcmc_arsvo(Y, X, Y0, dims, Hyper, sig2_hat, kappa1, kappa2, kappa3, kappa4, ...
    is_kappafixed, is_kappasym, nsim, burnin, pr)
% VAR_ARSVO_redu.m functionized line-for-line (clock-seed line 45 dropped).
% Draw order per sweep: alp -> beta -> h -> (mu,phi,sig2) -> kappa4, kappa1,
% kappa2 -> o -> po. Identical to mcmc_arsv except for the outlier scaling of
% the alp/beta/h blocks and the trailing o/po draw.
T = dims.T; n = dims.n; p = dims.p; k = dims.k; k_alp = dims.k_alp; k_beta = dims.k_beta;

ngrid = pr.svo.ngrid;                                   % [8]
o_grid = [1;linspace(pr.svo.o_grid_range(1),pr.svo.o_grid_range(2),ngrid)'];  % [9]

    % storage [12-19]
store_alp = zeros(nsim,k_alp);
store_beta = zeros(nsim,k_beta);
store_h = zeros(nsim,T,n);
store_o = zeros(nsim,T);
store_po = zeros(nsim,1);
store_hpara = zeros(nsim,3*n);
store_kappa = zeros(nsim,3);
count_phi = zeros(n,1);

    % initialize the Markov chain [22-42]
po = pr.svo.po_init;   % outlier probability
o = pr.svo.o_init*ones(T,1);
mu = log(sig2_hat);
sig2 = 1./gamrnd(Hyper.nuh,1./Hyper.Sh);
phi = Hyper.phi0;
A = zeros(k,n);
h = zeros(T,n);
XX = X'*X;
for ii=1:n
    iValpi = sparse(1:k,1:k,1./Hyper.Valp((ii-1)*k+1:ii*k));
    A(:,ii) = (XX + iValpi)\(X'*Y(:,ii));
    s2i = (Y(:,ii) - X*A(:,ii)).^2;
    mu(ii) = mean(log(s2i));
    h(:,ii) = bvar.sv.init_approx1N(s2i,mu(ii),phi(ii),sig2(ii));
end
alp = reshape(A,k_alp,1);                               %#ok<NASGU> % recomputed every sweep [37]
beta = zeros(k_beta,1);                                 %#ok<NASGU> % overwritten every sweep [38]
B0_id = nonzeros(tril(reshape(1:n^2,n,n),-1)');
B0 = eye(n);
[C_alp,idx_kappa1,idx_kappa2] = bvar.priors.minnesota_C(n,p,sig2_hat);
C_beta = Hyper.Vbeta/kappa4; % the constant part of Vbeta

for isim = 1:nsim + burnin
        % sample alp [49-65 -> bvar.samplers.eq_var_redu_tri with o]
    [~,Hyper.Valp] = bvar.priors.minn(p,kappa1,kappa2,kappa3,Y0,Y,4);
    U = zeros(T,n);                                     %#ok<NASGU> % dead in the legacy too [50]
    A = bvar.samplers.eq_var_redu_tri(Y,X,B0,h,A,Hyper.Valp,Hyper.alp0,o);
    alp = reshape(A,k_alp,1);

        % sample beta [68-81 -> bvar.samplers.alp_tri_cs with o; transposed]
    [~,Hyper.Vbeta] = bvar.priors.impact_B0(Y0,Y,kappa4);
    E = Y - X*A;
    beta = bvar.samplers.alp_tri_cs(E,h,Hyper.Vbeta,o)';
    B0(B0_id) = beta;

        % sample h [84-88]
    B0E = (E./o)*B0';
    for ii=1:n
        ystar = log(B0E(:,ii).^2 + pr.sv_offset);
        h(:,ii) = bvar.sv.ksc_ar1_mean(ystar,h(:,ii),mu(ii),phi(ii),sig2(ii));
    end

        % sample mu, phi and sig2 [91-92; ml_varsv phi bound .998]
    [mu,phi,sig2,flag_phi] = bvar.sv.sv_params(h,mu,phi,Hyper,pr.phi_mh_bnd);
    count_phi = count_phi + flag_phi;

        % sample kappa1, kappa2 and kappa4 [95-109]
    tmpc4 = sum(beta.^2./C_beta);
    kappa4 = gigrnd(Hyper.c0(3,1)-n*(n-1)/4,2*Hyper.c0(3,2),tmpc4,1);
    if is_kappafixed
        % nothing to do here
    elseif is_kappasym
        tmpc1 = sum(alp(idx_kappa1).^2./C_alp(idx_kappa1)) + ...
            sum(alp(idx_kappa2).^2./C_alp(idx_kappa2));
        kappa1 = gigrnd(Hyper.c0(1,1)-n^2*p/2,2*Hyper.c0(1,2),tmpc1,1);
        kappa2 = kappa1;
    else
        tmpc1 = sum(alp(idx_kappa1).^2./C_alp(idx_kappa1));
        tmpc2 = sum(alp(idx_kappa2).^2./C_alp(idx_kappa2));
        kappa1 = gigrnd(Hyper.c0(1,1)-n*p/2,2*Hyper.c0(1,2),tmpc1,1);
        kappa2 = gigrnd(Hyper.c0(2,1)-(n-1)*n*p/2,2*Hyper.c0(2,2),tmpc2,1);
    end

        % sample o and po [112-124 -> bvar.sv.svo_outlier]
    [o,po] = bvar.sv.svo_outlier(Y,X,A,B0,h,o_grid,po,Hyper.p0a,Hyper.p0b);

    if ( mod(isim, pr.progress_every) ==0 )             % [126-128]
        disp(  [ num2str(isim) ' loops... ' ] )
    end

    if isim > burnin                                    % [130-139]
        isave = isim-burnin;
        store_alp(isave,:) = reshape(A,1,k_alp);
        store_beta(isave,:) = beta';
        store_h(isave,:,:) = h;
        store_hpara(isave,:) = [mu',phi',sig2'];
        store_kappa(isave,:) = [kappa1,kappa2,kappa4];
        store_o(isave,:) = o;
        store_po(isave) = po;
    end
end

    % posterior summaries [143-150]
res = struct();
res.store_alp = store_alp;
res.store_beta = store_beta;
res.store_h = store_h;
res.store_o = store_o;
res.store_po = store_po;
res.store_hpara = store_hpara;
res.store_kappa = store_kappa;
res.count_phi = count_phi;
res.alp_mean = mean(store_alp)';
res.A_mean = reshape(res.alp_mean,k,n)';
res.beta_mean = mean(store_beta)';
res.h_mean = squeeze(mean(store_h));
res.hpara_mean = mean(store_hpara)';
res.kappa_mean = mean(store_kappa)';
res.o_mean = mean(store_o)';
res.po_mean = mean(store_po);
res.o_grid = o_grid;
end
