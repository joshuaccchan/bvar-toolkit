% bvar.forecast.iterate - one iterated-forecast + predictive-likelihood step for
% ONE posterior draw. Extracted 2026-09-01 (step 6, forecast engine). The entry
% point dispatches to internal NAMED branches whose bodies are VERBATIM copies
% of the legacy inline per-draw forecast blocks (same randn/rand/gamrnd call
% order and count, same expressions, same storage classes), so a caller that
% replays the legacy MCMC draw-for-draw and calls iterate once per kept draw
% reproduces the legacy tmpyhat arrays bitwise.
%
%   fc = bvar.forecast.iterate(branch, draw, cfg)
%
% branch - which legacy block to run (see the branch map below)
% draw   - struct of per-draw posterior quantities, LEGACY variable names
% cfg    - struct of vintage-level constants (fixed across draws)
% fc     - 2 x (2n+1) matrix, one row per legacy tmpyhat accumulator in the
%          branch, each row [EYtp1 (1 x n point forecast), lden (1 x n
%          per-variable log predictive likelihoods), lden_joint (joint log
%          predictive likelihood over ALL n modeled variables)]:
%            'mahp_sv'    : fc(1,:) -> tmpyhat1 (h=1), fc(2,:) -> tmpyhat4 (h=4)
%            'springer_*' : fc(1,:) -> tmpyhat0 (h=0), fc(2,:) -> tmpyhat1 (h=1)
%          A row whose evaluation guard is off this vintage (see "horizons"
%          below) comes back as zeros(1,2*n+1), matching the legacy convention
%          of leaving the preallocated zero row untouched; the simulation
%          draws for that step are still consumed, exactly as in the legacy
%          loop. The legacy expression sum(diag(log(CSig))) goes through a
%          complex intermediate when CSig has negative off-diagonal entries;
%          on R2025b diag() demotes the zero-imaginary diagonal back to real,
%          so real-data rows stay real (verified empirically), but on MATLAB
%          versions that retain the complex attribute the rows can come back
%          complex-typed with zero imaginary part. Either way the expression
%          is byte-identical to legacy on both sides - preserved deliberately
%          (do NOT "fix" it to log(diag(CSig))).
%
% Horizons: every branch simulates the legacy step loop in full (mahp_sv:
% tt = 1:4 evaluating at tt = 1 and 4; springer_*: tt = 1:2 evaluating at
% tt = 1 and 2, i.e. nowcast h=0 and one-quarter-ahead h=1). The first
% evaluated step is unguarded; each later step tt is evaluated only when
% t <= T - tt (verbatim legacy guard - note in the springer scripts this
% skips the h=1 evaluation at t = T-1 even though the outturn exists; that
% quirk is reproduced, not repaired).
%
% Inner-simulation count: all branches below use ONE simulated path per
% posterior draw (the legacy design - the predictive density is evaluated
% analytically conditional on the simulated volatility/lag path, and the
% outer MCMC loop provides the mixture). Families with genuine inner
% simulation loops (e.g. the OISV cluster forecasts) are NOT covered here
% and will add branches with an explicit count when canonicalized.
%
% Variable subsets: lden always covers the n variables OF THE MODEL and
% lden_joint all n jointly. Subsetting is the CALLER's job, done exactly as
% the legacy drivers do it: pass outturns already subsetted (springer model 1
% passes data_tpk(:,var_small) so the 4-variable model evaluates against the
% 4 outturn columns), and/or select columns downstream in the accumulation /
% table stage (bvar.forecast.tables).
%
% Missing-latest-observation convention (springer real-time vintages): when
% cfg.is_last_miss is true the branch first advances the state one extra
% simulation step (drawing volatility and Y) before the tt loop, so that
% "h=0" evaluates against the quarter after the last OBSERVED one - verbatim
% from the legacy scripts, including the rng draws consumed by the extra step.
%
% ---------------------------------------------------------------------------
% BRANCH MAP (which legacy inline blocks each branch canonicalizes; "body" =
% the per-draw forecast block copied verbatim, tmpyhat*(i,:) renamed fc(r,:)):
%
% 'mahp_sv'   chan2021_ijf_mahp/legacy/forecast_BVAR_MNG.m, body of the
%             forecast loop (lines 112-147; the CANONICAL copy). Textually
%             identical copies: forecast_BVAR_NG.m lines 113-148 and
%             forecast_BVAR_Minn.m lines 92-128. Structural BVAR with
%             per-variable random-walk SV: transforms (alp, beta, h_T, Sigh)
%             to reduced form, innovates all n log-volatilities each step.
%             draw: alp (k_alp x 1), beta (k_beta x 1), h_T (n x 1), Sigh (n x 1)
%             cfg:  Yt (estimation sample, last p rows feed the lag stack),
%                   Y (full outturn matrix; rows t+1 and t+4 are read),
%                   p, t, T
%
% 'springer_gauss'  chan2020_springer_largebvar/legacy/forecast_BVAR_Minn.m
%             lines 36-57 (CANONICAL). Textually identical bodies (given the
%             caller-supplied A/CSig/dSig and outturn subsetting):
%             forecast_BVAR_small.m 41-62 (caller passes data_tpk(:,var_small)),
%             forecast_BVAR_NCP.m 40-61, forecast_BVAR_IP.m 45-66,
%             forecast_BVAR_SSVS.m 52-73. Homoskedastic Gaussian errors.
%             draw: A (k x n; Minn/small/IP/SSVS callers reshape(beta,k,n)),
%                   CSig - IN THE LEGACY STORAGE CLASS: sparse(1:n,1:n,
%                   sqrt(Sig_hat)) for Minn/small, dense chol(Sig,'lower')
%                   for NCP/IP/SSVS (class changes sparse/full propagation),
%                   dSig (1 x n: Sig_hat' for Minn/small, diag(Sig)' for
%                   NCP/IP/SSVS - legacy recomputes it inside the step loop;
%                   it is constant there, so assigning the passed value in
%                   the same loop position is bit-identical)
%             cfg:  shortYt, data_tpk (>= 2 rows; pre-subsetted for model 1),
%                   is_last_miss, p, t, T
%
% 'springer_csv'  chan2020_springer_largebvar/legacy/forecast_BVAR_CSV.m
%             lines 67-94. Gaussian errors with common stochastic volatility
%             (CSV): scalar AR(1) log-volatility htp1 innovated before each
%             step, joint density carries -n/2*htp1 - .5*(u'u)/exp(htp1).
%             draw: A, CSig (dense chol(Sig,'lower')), Sig, h (Tt x 1 path;
%                   only h(end) is read, verbatim), rho, sigh2
%             cfg:  shortYt, data_tpk, is_last_miss, p, t, T
%
% 'springer_csv_t'  chan2020_springer_largebvar/legacy/forecast_BVAR_CSV_t.m
%             lines 77-106. CSV + Student-t errors: simulation divides by
%             sqrt(gamrnd(nu/2,2/nu)); densities are Student-t (ct/ct_joint).
%             draw: A, CSig, Sig, h, rho, sigh2, nu
%             cfg:  shortYt, data_tpk, is_last_miss, p, t, T
%
% 'springer_csv_t_ma'  chan2020_springer_largebvar/legacy/
%             forecast_BVAR_CSV_t_MA.m lines 109-142. CSV + t + MA(1) errors:
%             initializes the MA state from E = Hpsi\(shortYt - Z*A), adds
%             psi*etp1' to the conditional mean each step.
%             draw: A, CSig, Sig, h, rho, sigh2, nu, psi, Hpsi (the sparse
%                   Tt x Tt MA rotation held by the sampler)
%             cfg:  shortYt, Z (the Tt x k estimation design), data_tpk,
%                   is_last_miss, p, t, T
%
% NOT covered here (deferred to their family passes; the audit's remaining
% inline blocks): the 9 chan2020_jbes_kronecker/realtime_forecasts scripts
% (same skeleton at tt = 1:5 evaluating tt = 1,2,3,5, plus homoskedastic-t /
% MA / t-CSV / t-MA / CSV-MA / t-CSV-MA / n-variate-SV-small error families),
% cjz2019_ad_opthyper/forecast_BVAR_NCP.m (Gaussian, tt = 1:4 evaluating
% tt = 1 and 4, final-vintage data, no is_last_miss step), and the two
% chan_koop_yu2024_jbes_oisv cluster fragments forecast_CS_MH.m /
% forecast_SVARSV_MH.m (12-horizon monthly design, 2n+3-column rows).
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
% canonical: chan2021_ijf_mahp forecast_BVAR_MNG.m forecast-loop body 113-146
alp = draw.alp; beta = draw.beta; h_Tp1 = draw.h_T; Sigh = draw.Sigh;
Yt = cfg.Yt; Y = cfg.Y; p = cfg.p; t = cfg.t; T = cfg.T;
n = size(Y,2);
A_id = nonzeros(tril(reshape(1:n^2,n,n),-1)');  % forecast_BVAR_MNG.m line 11
A = eye(n);                                     % line 12; strict lower triangle fully overwritten below
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
% canonical: chan2020_springer_largebvar forecast_BVAR_Minn.m lines 36-57
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
    dSig = draw.dSig;   % legacy: dSig = Sig_hat' (Minn/small) / diag(Sig)' (NCP/IP/SSVS); constant across tt
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
% canonical: chan2020_springer_largebvar forecast_BVAR_CSV.m lines 67-94
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
% canonical: chan2020_springer_largebvar forecast_BVAR_CSV_t.m lines 77-106
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
% canonical: chan2020_springer_largebvar forecast_BVAR_CSV_t_MA.m lines 109-142
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
