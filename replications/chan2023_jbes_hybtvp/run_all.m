% chan2023_jbes_hybtvp/run_all - main_HYB_TVPSV.m, functionized: estimates the
% hybrid TVP-VAR of Chan (2023) and returns the stored draws together with the
% Savage-Dickey density ratios that produce the paper's log Bayes factors.
%
%   out = run_all()                                  % the legacy settings
%   out = run_all(nsim, burnin, seed, varid, is_gamfixed, is_kappafixed)
%
%   nsim, burnin  - default 50000 / 1000 (main_HYB_TVPSV.m 16-17)
%   seed          - seeds the stream ONCE, before initialization. The legacy file
%                   seeds from the clock at line 85, AFTER its initialize block
%                   has already drawn n volatility paths; removing that line (the
%                   only edit the equivalence test makes) leaves one continuous
%                   stream covering initialization and the loop, which is what
%                   this function reproduces.
%   varid         - variable selection; default preset.varid, the n = 6 set that
%                   is the active line in the legacy file. preset also carries
%                   varid_n3 and varid_n20.
%   is_gamfixed   - default false. When true, gam is held at [0,0; ones] rather
%                   than drawn (legacy 74-79).
%   is_kappafixed - default false (legacy 34).
%
% Output struct: the twelve store_* arrays under their legacy names, the data and
% design (Y, Y0, X2, varid, T, n, p, k_alp, k_beta), and the post-processing the
% legacy script prints - lBF (log Bayes factors against HYB(1,1), (1,0), (0,1),
% (0,0)), lpostgam, lprigam, gam_hat and gam_mode.
%
% Core used: bvar.priors.resid_var_allvars_ridge (legacy get_resid_var_v2),
% bvar.priors.minnesota_C (get_C), bvar.priors.vtheta (getVtheta - this package
% hard-codes kappa_3 = .2 and kappa_4 = 1, supplied here as kappa(3:4)),
% bvar.samplers.eq_hyb_tvp (sample_gam_thetai_ver2), bvar.sv.ksc_rw_h0
% (sample_SVRW - the two are bitwise identical, verified over 200 randomized
% inputs), bvar.util.gam_mode (get_gammode), third_party/gigrnd.m.
%
% Functionized 2026-09-05 (step 11). Draw-for-draw equivalence against the
% unmodified legacy script: tests/unit/test_hybtvp_equivalence.m.
%
% See:
% Chan, J.C.C. (2023). Large Hybrid Time-Varying Parameter VARs, Journal of
% Business and Economic Statistics, 41(3): 890-905.

function out = run_all(nsim, burnin, seed, varid, is_gamfixed, is_kappafixed)
thisdir = fileparts(mfilename('fullpath'));

    % make bvar.* and gigrnd resolvable when called standalone
if isempty(which('bvar.samplers.eq_hyb_tvp'))
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

if nargin < 1 || isempty(nsim),          nsim = pr.nsim;                   end
if nargin < 2 || isempty(burnin),        burnin = pr.burnin;               end
if nargin < 3 || isempty(seed),          seed = [];                        end
if nargin < 4 || isempty(varid),         varid = pr.varid;                 end
if nargin < 5 || isempty(is_gamfixed),   is_gamfixed = pr.is_gamfixed;     end
if nargin < 6 || isempty(is_kappafixed), is_kappafixed = pr.is_kappafixed; end
if ~isempty(seed)
    rng(seed, 'twister');
end

p = pr.p;

    % ---- data (main 21-31) ----
    % xlsread, not readmatrix, deliberately: on this file the two disagree in the
    % last bit of some cells (max absolute difference 6.9e-18 over the 238 x 248
    % block - they parse a few decimal strings to adjacent doubles). That is
    % numerically nothing, but it is enough to break the draw-for-draw guarantee
    % against the legacy script, which is the point of this driver. Keep the
    % legacy call even though xlsread is deprecated.
data = xlsread(fullfile(thisdir, 'legacy', pr.data_file), pr.data_range); %#ok<XLSRD>
Y0 = data(1:4, varid);      % first 4 obs as initial conditions
Y  = data(5:end, varid);
[T,n] = size(Y);
k_alp  = n*(n-1)/2;         % dimension of the impact matrix
k_beta = n^2*p + n;         % number of VAR coefficients

    % ---- prior (main 37-48) ----
kappa = pr.kappa0;
sig2 = bvar.priors.resid_var_allvars_ridge(Y0, Y);
[C, idx_kappa1, idx_kappa2] = bvar.priors.minnesota_C(n, p, sig2);
nuh0 = pr.nuh0_scalar*ones(n,1);
Sh0  = pr.Sh0_factor*(nuh0-1).*ones(n,1);
Vsigbeta = pr.Vsigbeta_offdiag*ones(k_beta,1);
Vsigbeta(1:n*p+1:end) = pr.Vsigbeta_intercept;       % intercepts
Vsigalp = pr.Vsigalp*ones(k_alp,1);
c01 = pr.c01;
c02 = pr.c02;
ah = pr.ah_scalar*ones(n,1);
Vh = pr.Vh_scalar*ones(n,1);
ap = pr.ap;
bp = pr.bp;

    % ---- design matrix (main 51-56) ----
tmpY = [Y0(end-p+1:end,:); Y];
X2 = zeros(T, n*p);
for ii = 1:p
    X2(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
end
X2 = [ones(T,1) X2];

    % ---- storage (main 59-71) ----
store_alp      = zeros(nsim,T,k_alp);
store_beta     = zeros(nsim,T,k_beta);
store_h        = zeros(nsim,T,n);
store_Sigbeta  = zeros(nsim,k_beta);
store_Sigalp   = zeros(nsim,k_alp);
store_Sigh     = zeros(nsim,n);
store_beta0    = zeros(nsim,k_beta);
store_alp0     = zeros(nsim,k_alp);
store_h0       = zeros(nsim,n);
store_p0       = zeros(nsim,n,2);
store_gam      = zeros(nsim,n*2);
store_kappa    = zeros(nsim,2);
store_lpostgam = zeros(nsim,4);

    % ---- initialize the chain (main 74-82 + initialize.m) ----
if is_gamfixed
    gam = [0,0; [zeros(n-1,1),ones(n-1,1)]];
else
    gam = zeros(n,2);
end
lp_gam = zeros(n,4);
lpostgam = zeros(1,4);
alp = zeros(T,k_alp);
beta = zeros(T,k_beta);

    % initialize.m, inline: beta0/alp0 by equation-wise least squares, then the
    % first volatility paths. This block DRAWS (n calls to the KSC sampler), and
    % in the legacy file it runs before the clock seeding on line 85.
beta0 = zeros(k_beta,1);
alp0  = zeros(k_alp,1);
count1 = 0;
U = zeros(T,n);
for ii = 1:n
    ki = n*p+1+ii-1;
    Xi = [X2 -Y(:,1:ii-1)];
    theta0 = (Xi'*Xi)\(Xi'*Y(:,ii));
    beta0((ii-1)*k_beta/n+1:ii*k_beta/n) = theta0(1:k_beta/n);
    alp0(count1+1:count1+ii-1) = theta0(k_beta/n+1:ki);
    U(:,ii) = Y(:,ii)-Xi*theta0;
    count1 = count1 + ii-1;
end
Sigalp  = Vsigalp;
Sigbeta = Vsigbeta;
p0 = pr.p0_init*ones(n,2);
h0 = mean(log(U.^2))';
Sigh = Sh0;
h = repmat(h0',T,1);
for ii = 1:n
    h(:,ii) = bvar.sv.ksc_rw_h0(log(U(:,ii).^2 + pr.sv_offset_init), h(:,ii), Sigh(ii), h0(ii));
end
    % U now holds the equation-wise least-squares residuals, as it does in the
    % legacy file: main 80 zeroes it, initialize.m then fills it. The loop only
    % ever writes U, so its value here does not enter the chain.

    % ---- MCMC (main 89-184) ----
disp('Starting MCMC for the hybrid TVP-VAR...');
start_time = clock; %#ok<CLOCK> % legacy timing call, kept

for isim = 1:nsim + burnin
        % sample states and parameters equation by equation
    [Valp,Vbeta] = bvar.priors.vtheta(idx_kappa1, idx_kappa2, ...
        [kappa(1), kappa(2), pr.kappa_3, pr.kappa_4], C, sig2);
    for ii = 1:n
        ki = n*p+1+ii-1;
        idx_a0 = (ii-1)*(ii-2)/2+1; idx_a1 = ii*(ii-1)/2;  % index for alp
        idx_b0 = (ii-1)*k_beta/n+1; idx_b1 = ii*k_beta/n;  % index for beta
        Xi = [X2 -Y(:,1:ii-1)];
        Yi = Y(:,ii);
        hi = h(:,ii);
        gami = gam(ii,:);
        p0i = p0(ii,:)';
        SigBeta = reshape(Sigbeta,k_beta/n,n);
        Sigbetai = SigBeta(:,ii);
        Sigalpi = Sigalp(idx_a0:idx_a1);

            % sample gami and thetai
        Sigthetai = [Sigbetai; Sigalpi];
        thetai0 = [beta0(idx_b0:idx_b1); alp0(idx_a0:idx_a1)];
        [gami,thetai,Ui,tilde_thetai,lp_gami] = bvar.samplers.eq_hyb_tvp(Xi,...
            thetai0,Sigthetai,Yi,hi,p0i,n*p+1,is_gamfixed,gami);
        gam(ii,:) = gami;
        lp_gam(ii,:) = lp_gami;
        Thetai = reshape(thetai,ki,T)';
        beta(:,idx_b0:idx_b1) = Thetai(:,1:k_beta/n);
        alp(:,idx_a0:idx_a1) = Thetai(:,k_beta/n+1:end);
        U(:,ii) = Ui;

            % sample hi  (note: no 1e-4 offset here, unlike initialization)
        hi = bvar.sv.ksc_rw_h0(log(Ui.^2),hi,Sigh(ii),h0(ii));
        h(:,ii) = hi;

            % sample thetai0 and Sigthetai
        tilde_Thetai = reshape(tilde_thetai,ki,T)';
        Vmui = [Vbeta(idx_b0:idx_b1); Valp(idx_a0:idx_a1); Vsigbeta(idx_b0:idx_b1); Vsigalp(idx_a0:idx_a1)];
        Wi = [Xi Xi.*tilde_Thetai];
        WiSig  = Wi'*sparse(1:T,1:T,exp(-hi));
        Kmui = sparse(1:2*ki,1:2*ki,1./Vmui) + WiSig*Wi;
        mui_hat = Kmui\(WiSig*Yi);
        mui = mui_hat + chol(Kmui,'lower')'\randn(2*ki,1);
        beta0(idx_b0:idx_b1) = mui(1:k_beta/n);
        alp0(idx_a0:idx_a1) = mui(k_beta/n+1:ki);
        Sigbeta(idx_b0:idx_b1) = mui(ki+1:ki+k_beta/n).^2;
        Sigalp(idx_a0:idx_a1) = mui(k_beta/n+ki+1:end).^2;
    end

        % sample h0
    Kh0 = sparse(1:n,1:n,1./Sigh + 1./Vh);
    h0_hat = Kh0\(ah./Vh + h(1,:)'./Sigh);
    h0 = h0_hat + chol(Kh0,'lower')'\randn(n,1);

        % sample Sigh
    eh = h - [h0';h(1:T-1,:)];
    Sigh = 1./gamrnd(nuh0+T/2,1./(Sh0 + sum(eh.^2)'/2));

        % sample kappa1 and kappa2
    if ~is_kappafixed
        tmpc1 = sum(beta0(idx_kappa1).^2./C(idx_kappa1));
        tmpc2 = sum(beta0(idx_kappa2).^2./C(idx_kappa2));
        kappa(1) = gigrnd(c01(1)-n*p/2,2*c01(2),tmpc1,1);
        kappa(2) = gigrnd(c02(1)-(n-1)*n*p/2,2*c02(2),tmpc2,1);
    end

        % sample p0
    if ~is_gamfixed
        p0 = betarnd(ap+gam,bp+1-gam);
    end

        % evaluate posterior of gamma
    lpostgam(1) = sum(lp_gam(:,1));                   % HYB-(0,0)
    lpostgam(2) = lp_gam(1,1) + sum(lp_gam(2:end,2)); % HYB-(0,1)
    lpostgam(3) = sum(lp_gam(:,3));                   % HYB-(1,0)
    lpostgam(4) = lp_gam(1,3)+sum(lp_gam(2:end,4));   % HYB-(1,1)

    if isim > burnin
        isave = isim - burnin;
        store_h(isave,:,:) = h;
        store_alp(isave,:,:) = alp;
        store_beta(isave,:,:) = beta;
        store_Sigbeta(isave,:) = Sigbeta';
        store_Sigalp(isave,:) = Sigalp';
        store_Sigh(isave,:) = Sigh';
        store_beta0(isave,:) = beta0';
        store_alp0(isave,:) = alp0';
        store_h0(isave,:) = h0';
        store_p0(isave,:,:) = p0;
        store_gam(isim-burnin,:) = reshape(gam',1,n*2);
        store_kappa(isave,:) = kappa';
        store_lpostgam(isave,:) = lpostgam;
    end

    if (mod(isim, pr.progress_every) == 0)
        disp([num2str(isim) ' loops... ']);
    end
end

disp( ['MCMC takes '  num2str( etime( clock, start_time) ) ' seconds' ] ); %#ok<CLOCK,DETIM>
disp(' ' );

    % ---- SDDR / log Bayes factors (main 190-200) ----
maxtmp = max(store_lpostgam);
lpostgam_hat = log(mean(exp(store_lpostgam-repmat(maxtmp,nsim,1)))) + maxtmp;
lprigam_c = n*sum(gammaln(ap+bp)-gammaln(ap+bp+1)-gammaln(ap)-gammaln(bp));
lprigam = zeros(1,4);
lprigam(1) = lprigam_c + n*sum(gammaln(ap)+gammaln(bp+1));
lprigam(2) = lprigam_c + n*(gammaln(ap(1)+1)+gammaln(bp(1))+gammaln(ap(2))+gammaln(bp(2)+1));
lprigam(3) = lprigam_c + n*(gammaln(ap(1))+gammaln(bp(1)+1)+gammaln(ap(2)+1)+gammaln(bp(2)));
lprigam(4) = lprigam_c + n*sum(gammaln(ap+1)+gammaln(bp));
lBF = lprigam - lpostgam_hat;
gam_mode = reshape(bvar.util.gam_mode(store_gam),2,n)';
gam_hat = reshape(mean(store_gam),2,n)';

out = struct();
out.Y = Y; out.Y0 = Y0; out.X2 = X2; out.varid = varid;
out.T = T; out.n = n; out.p = p; out.k_alp = k_alp; out.k_beta = k_beta;
out.nsim = nsim; out.burnin = burnin; out.seed = seed;
out.is_gamfixed = is_gamfixed; out.is_kappafixed = is_kappafixed;
out.store_alp = store_alp;         out.store_beta = store_beta;
out.store_h = store_h;             out.store_Sigbeta = store_Sigbeta;
out.store_Sigalp = store_Sigalp;   out.store_Sigh = store_Sigh;
out.store_beta0 = store_beta0;     out.store_alp0 = store_alp0;
out.store_h0 = store_h0;           out.store_p0 = store_p0;
out.store_gam = store_gam;         out.store_kappa = store_kappa;
out.store_lpostgam = store_lpostgam;
out.lpostgam = lpostgam_hat; out.lprigam = lprigam; out.lBF = lBF;
out.gam_hat = gam_hat; out.gam_mode = gam_mode;
out.sig2 = sig2; out.C = C; out.idx_kappa1 = idx_kappa1; out.idx_kappa2 = idx_kappa2;
out.preset = pr;
end
