% bvar.forecast.iterate - one iterated-forecast + predictive-likelihood step for
% ONE posterior draw. The entry point dispatches to internal NAMED branches,
% one per error/volatility specification; call it once per kept draw.
%
%   fc = bvar.forecast.iterate(branch, draw, cfg)
%
% branch - which specification to run (see the branch list below)
% draw   - struct of per-draw posterior quantities; the fields each branch
%          reads are listed with that branch
% cfg    - struct of vintage-level constants (fixed across draws)
% fc     - 2 x (2n+1) matrix, one row per horizon the branch evaluates, each
%          row [EYtp1 (1 x n point forecast), lden (1 x n per-variable log
%          predictive likelihoods), lden_joint (joint log predictive
%          likelihood over ALL n modeled variables)]:
%            'mahp_sv'    : fc(1,:) -> h=1, fc(2,:) -> h=4
%            'springer_*' : fc(1,:) -> h=0 (nowcast), fc(2,:) -> h=1
%          A row whose evaluation guard is off this vintage (see "Horizons"
%          below) comes back as zeros(1,2*n+1); the simulation draws for that
%          step are still consumed.
%
% TRAP (complex-typed rows): the expression sum(diag(log(CSig))) goes through
% a complex intermediate when CSig has negative off-diagonal entries. On
% R2025b diag() demotes the zero-imaginary diagonal back to real, so real-data
% rows stay real (verified empirically); on MATLAB versions that retain the
% complex attribute the rows can come back complex-typed with zero imaginary
% part, and a downstream max() then compares by magnitude. The expression is
% preserved deliberately - do NOT "fix" it to log(diag(CSig)).
%
% Horizons: every branch simulates the step loop in full (mahp_sv: tt = 1:4,
% evaluating at tt = 1 and 4; springer_*: tt = 1:2, evaluating at tt = 1 and
% 2, i.e. nowcast h=0 and one-quarter-ahead h=1). The first evaluated step is
% unguarded; each later step tt is evaluated only when t <= T - tt. In the
% springer branches that guard skips the h=1 evaluation at t = T-1 even though
% the outturn exists; the quirk is preserved deliberately.
%
% Inner-simulation count: all branches below use ONE simulated path per
% posterior draw - the predictive density is evaluated analytically
% conditional on the simulated volatility/lag path, and the outer MCMC loop
% provides the mixture. Families with genuine inner simulation loops (e.g. the
% OISV cluster forecasts) are NOT covered here and will add branches with an
% explicit count when canonicalized.
%
% Variable subsets: lden always covers the n variables OF THE MODEL and
% lden_joint all n jointly. Subsetting is the CALLER's job: pass outturns
% already subsetted (a 4-variable model passes data_tpk(:,var_small), so it
% evaluates against the 4 outturn columns), and/or select columns downstream
% in the accumulation / table stage (bvar.forecast.tables).
%
% Missing-latest-observation convention (real-time vintages): when
% cfg.is_last_miss is true the branch first advances the state one extra
% simulation step (drawing volatility and Y) before the tt loop, so that "h=0"
% evaluates against the quarter after the last OBSERVED one. That extra step
% consumes rng draws.
%
% ---------------------------------------------------------------------------
% BRANCHES: what each one models, and the draw/cfg fields it reads.
%
% 'mahp_sv'   Structural BVAR with per-variable random-walk SV: transforms
%             (alp, beta, h_T, Sigh) to reduced form, innovates all n
%             log-volatilities each step.
%             draw: alp (k_alp x 1), beta (k_beta x 1), h_T (n x 1), Sigh (n x 1)
%             cfg:  Yt (estimation sample, last p rows feed the lag stack),
%                   Y (full outturn matrix; rows t+1 and t+4 are read),
%                   p, t, T
%
% 'springer_gauss'  Homoskedastic Gaussian errors.
%             draw: A (k x n coefficient matrix; a caller holding
%                     beta = vec(A) passes reshape(beta,k,n))
%                   CSig - the error factor, IN THE STORAGE CLASS THE MODEL
%                     USES: sparse(1:n,1:n,sqrt(Sig_hat)) when Sig is
%                     diagonal, dense chol(Sig,'lower') otherwise. The class
%                     changes sparse/full propagation downstream, so pass the
%                     one the model actually uses
%                   dSig (1 x n) - the diagonal of Sig as a row: Sig_hat' in
%                     the diagonal case, diag(Sig)' otherwise. Held constant
%                     across the step loop
%             cfg:  shortYt, data_tpk (>= 2 rows; pre-subsetted when the model
%                   covers a subset of the outturn columns),
%                   is_last_miss, p, t, T
%
% 'springer_csv'  Gaussian errors with common stochastic volatility (CSV):
%             scalar AR(1) log-volatility htp1 innovated before each step,
%             joint density carries -n/2*htp1 - .5*(u'u)/exp(htp1).
%             draw: A, CSig (dense chol(Sig,'lower')), Sig, h (Tt x 1 path;
%                   only h(end) is read), rho, sigh2
%             cfg:  shortYt, data_tpk, is_last_miss, p, t, T
%
% 'springer_csv_t'  CSV + Student-t errors: simulation divides by
%             sqrt(gamrnd(nu/2,2/nu)); densities are Student-t (ct/ct_joint).
%             draw: A, CSig, Sig, h, rho, sigh2, nu
%             cfg:  shortYt, data_tpk, is_last_miss, p, t, T
%
% 'springer_csv_t_ma'  CSV + t + MA(1) errors: initializes the MA state from
%             E = Hpsi\(shortYt - Z*A), adds psi*etp1' to the conditional mean
%             each step.
%             draw: A, CSig, Sig, h, rho, sigh2, nu, psi, Hpsi (the sparse
%                   Tt x Tt MA rotation held by the sampler)
%             cfg:  shortYt, Z (the Tt x k estimation design), data_tpk,
%                   is_last_miss, p, t, T
% ---------------------------------------------------------------------------
%
% See:
% Chan, J.C.C. (2021). Minnesota-Type Adaptive Hierarchical Priors for
% Large Bayesian VARs, International Journal of Forecasting, 37(3), 1212-1226.
% Chan, J.C.C. (2020). Large Bayesian Vector Autoregressions. In: P. Fuleky (Ed.),
% Macroeconomic Forecasting in the Era of Big Data, 95-125, Springer, Cham.

function fc = iterate(branch, draw, cfg)
switch branch
    case 'mahp_sv'
        fc = mahp_sv(draw, cfg);
    case 'springer_gauss'
        fc = springer_gauss(draw, cfg);
    case 'springer_csv'
        fc = springer_csv(draw, cfg);
    case 'springer_csv_t'
        fc = springer_csv_t(draw, cfg);
    case 'springer_csv_t_ma'
        fc = springer_csv_t_ma(draw, cfg);
    otherwise
        error('bvar:forecast:iterate:unknownBranch', ...
            ['unknown branch ''%s''; use mahp_sv, springer_gauss, ' ...
            'springer_csv, springer_csv_t or springer_csv_t_ma'], branch);
end
end

% ---------------------------------------------------------------------------
function fc = mahp_sv(draw, cfg)
alp = draw.alp; beta = draw.beta; h_Tp1 = draw.h_T; Sigh = draw.Sigh;
Yt = cfg.Yt; Y = cfg.Y; p = cfg.p; t = cfg.t; T = cfg.T;
n = size(Y,2);
A_id = nonzeros(tril(reshape(1:n^2,n,n),-1)');  % strict lower triangle indices
A = eye(n);                                     % strict lower triangle fully overwritten below
fc = zeros(2,2*n+1);
    % trasnform the parameters into reduced-form
sqrtSigh = sqrt(Sigh);
A(A_id) = alp;
h_Tp1 = h_Tp1 + sqrtSigh.*randn(n,1);
S = (A\sparse(1:n,1:n,exp(h_Tp1)))/A';
B = (A\(reshape(beta,n*p+1,n)'))';
xtp1 = [1 reshape(Yt(end:-1:end-p+1,:)',1,n*p)];
CS = chol(S,'lower');
for tt=1:4
    EYtp1 = xtp1*B;
    dS = diag(S)';
    if tt == 1
        tmpu = CS\(Y(t+1,:)-EYtp1)';
        lden_joint = -n/2*log(2*pi) -sum(diag(log(CS))) -.5*(tmpu'*tmpu);
        lden = -.5*log(2*pi*dS) - .5*(Y(t+1,:)-EYtp1).^2./dS;
        fc(1,:) = [EYtp1 lden lden_joint];
    elseif tt == 4 && t<=T-tt
        tmpu = CS\(Y(t+4,:)-EYtp1)';
        lden_joint = -n/2*log(2*pi) -sum(diag(log(CS))) -.5*(tmpu'*tmpu);
        lden = -.5*log(2*pi*dS) - .5*(Y(t+4,:)-EYtp1).^2./dS;
        fc(2,:) = [EYtp1 lden lden_joint];
    end
    Ytp1 = EYtp1 + (CS*randn(n,1))';
    xtp1 = [1 Ytp1 xtp1(2:end-n)];

    h_Tp1 = h_Tp1 + sqrtSigh.*randn(n,1);
    S = (A\sparse(1:n,1:n,exp(h_Tp1)))/A';
    CS = chol(S,'lower');
end
end

% ---------------------------------------------------------------------------
function fc = springer_gauss(draw, cfg)
A = draw.A; CSig = draw.CSig;
shortYt = cfg.shortYt; data_tpk = cfg.data_tpk; is_last_miss = cfg.is_last_miss;
p = cfg.p; t = cfg.t; T = cfg.T;
n = size(A,2);
fc = zeros(2,2*n+1);
xtp1 = [1 reshape(shortYt(end:-1:end-p+1,:)',1,n*p)];
if is_last_miss % if the lastest data are missing, do one more forecast horizon
    Ytp1 = xtp1*A + (CSig*randn(n,1))';
    xtp1 = [1 Ytp1 xtp1(2:end-n)];
end
for tt=1:2
    EYtp1 = xtp1*A;
    dSig = draw.dSig;   % the diagonal of Sig as a row; constant across tt
    if tt == 1
        tmpu = CSig\(data_tpk(1,:)-EYtp1)';
        lden_joint = -n/2*log(2*pi) -sum(diag(log(CSig))) -.5*(tmpu'*tmpu);
        lden = -.5*log(2*pi*dSig) - .5*(data_tpk(1,:)-EYtp1).^2./dSig;
        fc(1,:) = [EYtp1 lden lden_joint];
    elseif tt == 2 && t<=T-tt
        tmpu = CSig\(data_tpk(2,:)-EYtp1)';
        lden_joint = -n/2*log(2*pi) -sum(diag(log(CSig))) -.5*(tmpu'*tmpu);
        lden = -.5*log(2*pi*dSig) - .5*(data_tpk(2,:)-EYtp1).^2./dSig;
        fc(2,:) = [EYtp1 lden lden_joint];
    end
    Ytp1 = EYtp1 + (CSig*randn(n,1))';
    xtp1 = [1 Ytp1 xtp1(2:end-n)];
end
end

% ---------------------------------------------------------------------------
function fc = springer_csv(draw, cfg)
A = draw.A; CSig = draw.CSig; Sig = draw.Sig; h = draw.h;
rho = draw.rho; sigh2 = draw.sigh2;
shortYt = cfg.shortYt; data_tpk = cfg.data_tpk; is_last_miss = cfg.is_last_miss;
p = cfg.p; t = cfg.t; T = cfg.T;
n = size(A,2);
fc = zeros(2,2*n+1);
xtp1 = [1 reshape(shortYt(end:-1:end-p+1,:)',1,n*p)];
htp1 = h(end);
if is_last_miss % if the lastest data are missing, do one more forecast horizon
    htp1 = rho*htp1 + sqrt(sigh2)*randn;
    EYtp1 = xtp1*A;
    Ytp1 = EYtp1 + (exp(htp1/2)*CSig*randn(n,1))';
    xtp1 = [1 Ytp1 xtp1(2:end-n)];
end
for tt = 1:2
    htp1 = rho*htp1 + sqrt(sigh2)*randn;
    EYtp1 = xtp1*A;
    dSig = exp(htp1)*diag(Sig)';
    if tt == 1
        tmpu = CSig\(data_tpk(1,:)-EYtp1)';
        lden_joint = -n/2*log(2*pi) -n/2*htp1 -sum(diag(log(CSig)))...
            -.5*(tmpu'*tmpu)/exp(htp1);
        lden = -.5*log(2*pi*dSig) - .5*(data_tpk(1,:)-EYtp1).^2./dSig;
        fc(1,:) = [EYtp1 lden lden_joint];
    elseif tt == 2 && t<=T-tt
        tmpu = CSig\(data_tpk(2,:)-EYtp1)';
        lden_joint = -n/2*log(2*pi) -n/2*htp1 -sum(diag(log(CSig)))...
            -.5*(tmpu'*tmpu)/exp(htp1);
        lden = -.5*log(2*pi*dSig) - .5*(data_tpk(2,:)-EYtp1).^2./dSig;
        fc(2,:) = [EYtp1 lden lden_joint];
    end
    Ytp1 = EYtp1 + (exp(htp1/2)*CSig*randn(n,1))';
    xtp1 = [1 Ytp1 xtp1(2:end-n)];
end
end

% ---------------------------------------------------------------------------
function fc = springer_csv_t(draw, cfg)
A = draw.A; CSig = draw.CSig; Sig = draw.Sig; h = draw.h;
rho = draw.rho; sigh2 = draw.sigh2; nu = draw.nu;
shortYt = cfg.shortYt; data_tpk = cfg.data_tpk; is_last_miss = cfg.is_last_miss;
p = cfg.p; t = cfg.t; T = cfg.T;
n = size(A,2);
fc = zeros(2,2*n+1);
xtp1 = [1 reshape(shortYt(end:-1:end-p+1,:)',1,n*p)];
htp1 = h(end);
if is_last_miss % if the lastest data are missing, do one more forecast horizon
    htp1 = rho*htp1 + sqrt(sigh2)*randn;
    EYtp1 = xtp1*A;
    Ytp1 = EYtp1 + (exp(htp1/2)*CSig*randn(n,1))'/sqrt(gamrnd(nu/2,2/nu));
    xtp1 = [1 Ytp1 xtp1(2:end-n)];
end
for tt = 1:2
    htp1 = rho*htp1 + sqrt(sigh2)*randn;
    EYtp1 = xtp1*A;
    dSig = exp(htp1)*diag(Sig)';
    ct = gammaln((nu+1)/2) - gammaln(nu/2) - .5*log(nu*pi*dSig);
    ct_joint = gammaln((nu+n)/2) - gammaln(nu/2) - n/2*log(nu*pi);
    if tt == 1
        tmpu = CSig\(data_tpk(1,:)-EYtp1)';
        lden_joint = ct_joint - sum(diag(log(CSig))) -n/2*htp1...
            -(nu+n)/2*log(1 + (tmpu'*tmpu)/nu/exp(htp1));
        lden = ct - (nu+1)/2*log(1 + (data_tpk(1,:)-EYtp1).^2./dSig/nu);
        fc(1,:) = [EYtp1 lden lden_joint];
    elseif tt == 2 && t<=T-tt
        tmpu = CSig\(data_tpk(2,:)-EYtp1)';
        lden_joint = ct_joint - sum(diag(log(CSig))) -n/2*htp1...
            -(nu+n)/2*log(1 + (tmpu'*tmpu)/nu/exp(htp1));
        lden = ct - (nu+1)/2*log(1 + (data_tpk(2,:)-EYtp1).^2./dSig/nu);
        fc(2,:) = [EYtp1 lden lden_joint];
    end
    Ytp1 = EYtp1 + (exp(htp1/2)*CSig*randn(n,1))'/sqrt(gamrnd(nu/2,2/nu));
    xtp1 = [1 Ytp1 xtp1(2:end-n)];
end
end

% ---------------------------------------------------------------------------
function fc = springer_csv_t_ma(draw, cfg)
A = draw.A; CSig = draw.CSig; Sig = draw.Sig; h = draw.h;
rho = draw.rho; sigh2 = draw.sigh2; nu = draw.nu; psi = draw.psi; Hpsi = draw.Hpsi;
shortYt = cfg.shortYt; Z = cfg.Z; data_tpk = cfg.data_tpk; is_last_miss = cfg.is_last_miss;
p = cfg.p; t = cfg.t; T = cfg.T;
n = size(A,2);
fc = zeros(2,2*n+1);
xtp1 = [1 reshape(shortYt(end:-1:end-p+1,:)',1,n*p)];
E = Hpsi\(shortYt - Z*A);
etp1 = E(end,:)';
htp1 = h(end);
if is_last_miss % if the lastest data are missing, do one more forecast horizon
    htp1 = rho*htp1 + sqrt(sigh2)*randn;
    EYtp1 = xtp1*A + psi*etp1';
    etp1 = (exp(htp1/2)*CSig*randn(n,1))/sqrt(gamrnd(nu/2,2/nu));
    Ytp1 = EYtp1 + etp1';
    xtp1 = [1 Ytp1 xtp1(2:end-n)];
end
for tt = 1:2
    htp1 = rho*htp1 + sqrt(sigh2)*randn;
    EYtp1 = xtp1*A + psi*etp1';
    dSig = exp(htp1)*diag(Sig)';
    ct = gammaln((nu+1)/2) - gammaln(nu/2) - .5*log(nu*pi*dSig);
    ct_joint = gammaln((nu+n)/2) - gammaln(nu/2) - n/2*log(nu*pi);
    if tt == 1
        tmpu = CSig\(data_tpk(1,:)-EYtp1)';
        lden_joint = ct_joint - sum(diag(log(CSig))) -n/2*htp1...
            -(nu+n)/2*log(1 + (tmpu'*tmpu)/nu/exp(htp1));
        lden = ct - (nu+1)/2*log(1 + (data_tpk(1,:)-EYtp1).^2./dSig/nu);
        fc(1,:) = [EYtp1 lden lden_joint];
    elseif tt == 2 && t<=T-tt
        tmpu = CSig\(data_tpk(2,:)-EYtp1)';
        lden_joint = ct_joint - sum(diag(log(CSig))) -n/2*htp1...
            -(nu+n)/2*log(1 + (tmpu'*tmpu)/nu/exp(htp1));
        lden = ct - (nu+1)/2*log(1 + (data_tpk(2,:)-EYtp1).^2./dSig/nu);
        fc(2,:) = [EYtp1 lden lden_joint];
    end
    etp1 = (exp(htp1/2)*CSig*randn(n,1))/sqrt(gamrnd(nu/2,2/nu));
    Ytp1 = EYtp1 + etp1';
    xtp1 = [1 Ytp1 xtp1(2:end-n)];
end
end
