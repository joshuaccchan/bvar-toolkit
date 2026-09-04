% bvar.ml.kron_bvar_t_ma - log marginal likelihood of the BVAR-t-MA model of
% Chan (2020, JBES): analytic Student-t MA(1) likelihood at the posterior
% means (lam integrated analytically), prior ordinate, posterior ordinates
% for (A,Sig) and nu Rao-Blackwellized over the stored (lam, psi) draws, and
% a psi ordinate from a reduced run of nsims sweeps re-drawing (lam, psi) at
% fixed (A_mean, Sig_mean, nu_mean), continuing from the final stored psi
% draw with the psi-MH proposal warm-started at the ESTIMATION run's final
% mode psihat. Consumes rng in the reduced run only (gamrnd, randn, rand);
% fminunc/fminbnd are deterministic.
%
% Extracted 2026-09-02 from chan2020_jbes_kronecker/legacy/ml_BVAR_t_MA.m,
% body verbatim. Clean bill: every ordinate sits at the same starred point -
% in particular line 19 uses s2(1)/(1+psi_mean^2), the term ml_BVAR_MA.m
% gets wrong. Leftover-workspace reads are chain/optimizer continuation made
% explicit (store_theta(nsims,1), est.state.psihat, reconstructed optimset).
% No bugcompat flag needed; verbatim quirks are in tests/variant_map.md.
%
%   [ML, out] = bvar.ml.kron_bvar_t_ma(shortY, X, pri, est)
%
%   pri: A0, VA0, nu0, S0, psi0, Vpsi, nuub
%   est: nsims, store_A, store_Sig (running sums), store_lam,
%        store_theta ([psi nu] columns), state.psihat
%   out: llike, lpri, lpost, store_lpost (reduced-run den_psi column),
%        store_lpost1, A_mean, Sig_mean, theta_mean
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics, 38(1), 68-79.

function [ML, out] = kron_bvar_t_ma(shortY, X, pri, est)
[T, n] = size(shortY);
k = size(X, 2);
A0 = pri.A0; VA0 = pri.VA0; nu0 = pri.nu0; S0 = pri.S0;
psi0 = pri.psi0; Vpsi = pri.Vpsi; nuub = pri.nuub;
nsims = est.nsims;
store_A = est.store_A; store_Sig = est.store_Sig;
store_lam = est.store_lam; store_theta = est.store_theta;
psihat = est.state.psihat;              % estimation run's final psi-MH mode
options = optimset('Display', 'off', 'LargeScale','off') ;  % = BVAR_t_MA.m line 28
lpri_psi = @(x) -.5*(x-psi0)^2/Vpsi -1e10*(x<-.99 || x>.99);    % = BVAR_t_MA.m line 10

    % estimation-tail posterior means [BVAR_t_MA.m lines 137-139]
A_mean = store_A/nsims;
Sig_mean = store_Sig/nsims;
theta_mean = mean(store_theta)';

    % [ml_BVAR_t_MA.m lines 9-27]
psi_mean = theta_mean(1);
nu_mean = theta_mean(2);
ngrid = 300;

    % evaluate the log likelihood
Hpsi = speye(T) + psi_mean(1)*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
CSig = chol(Sig_mean,'lower');
Utld = Hpsi\(shortY - X*A_mean);
tmp = (Utld/CSig');
s2 = sum(tmp.^2,2);
s2(1) = s2(1)/(1+psi_mean^2);
llike = -T*n/2*log(nu_mean*pi) - n/2*log(1+psi_mean^2) ...
    + T*(gammaln((nu_mean+n)/2) - gammaln(nu_mean/2)) ...
    - T*sum(log(diag(CSig))) - (nu_mean+n)/2*sum(log(1+s2/nu_mean));

c_psi = 1/(normcdf(1,psi0,sqrt(Vpsi))-normcdf(-1,psi0,sqrt(Vpsi)));
lpri = bvar.ml.lniwpdf(A_mean,Sig_mean,A0,sparse(1:k,1:k,1./VA0),nu0,S0) ...
    + log(1/(nuub-2)) ...
    -.5*log(2*pi*Vpsi) + log(c_psi) -.5*(psi_mean(1)-psi0)^2/Vpsi;

    % evaluate the posterior density [lines 30-64]
store_lpost = zeros(nsims,2); % [log density of A Sig, density of nu]
nugrid = sort([nu_mean; linspace(2,nuub,ngrid)']);
nuidx = find(nugrid==nu_mean);
for isim = 1:nsims
    lam = store_lam(isim,:)';
    psi1 = store_theta(isim,1);
    Hpsi = speye(T) + psi1*sparse(2:T,1:(T-1),ones(1,T-1),T,T);

        % compute the conditional density of Sig and A
    Xtld = Hpsi\X;
    Ytld = Hpsi\shortY;
    iO_lam = sparse(1:T,1:T,1./lam);
    XiO = Xtld'*iO_lam;
    KA = sparse(1:k,1:k,1./VA0) + XiO*Xtld;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiO*Ytld);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + Ytld'*iO_lam*Ytld ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    lden_ASig = bvar.ml.lniwpdf(A_mean,Sig_mean,Ahat,KA,nu0+T,Shat);

        % compute the conditional density of nu
    sum1 = sum(log(lam));
    sum2 = sum(1./lam);
    fnu = @(x) T*(x/2.*log(x/2)-gammaln(x/2)) - (x/2+1)*sum1 - x/2*sum2;
    tmpden = fnu(nugrid);
    tmpden = exp(tmpden-max(tmpden));
    tmpden = tmpden/(sum(tmpden)*(nugrid(2)-nugrid(1)));
    den_nu = tmpden(nuidx);

    store_lpost(isim,:) = [lden_ASig den_nu];
end
tmpmax = max(store_lpost(:,1));
lpost = zeros(3,1);
lpost(1) = log(mean(exp(store_lpost(:,1)-tmpmax))) + tmpmax;
lpost(2) = log(mean(store_lpost(:,2)));

    % psi ordinate: reduced run [lines 66-118]; psi1/Hpsi continue from the
    % last stored draw, psihat from the estimation run's final mode
store_lpost1 = store_lpost;
store_lpost = zeros(nsims,1); % density of psi
psigrid = sort([theta_mean(1); linspace(-.99,.99,ngrid)']);
psiidx = find(psigrid==theta_mean(1));

U = shortY - X*A_mean;
CSig = chol(Sig_mean,'lower');
for isim = 1:nsims

    % sample lam
    Utld = Hpsi\U;
    tmp = Utld/CSig';
    s2 = sum(tmp.^2,2);
    s2(1) = s2(1)/(1+psi1^2);
    lam = 1./gamrnd((n+nu_mean)/2,2./(s2+nu_mean));

    % sample psi
    lp_psi = @(x) bvar.ml.llike_csv_ma(x,U,Sig_mean,log(lam)) + lpri_psi(x);
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
        Hpsi = speye(T) + psi1*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
    end

        % compute the conditional density of psi [normalization centered at
        % tmpden(psiidx) - verbatim m6 quirk, mathematically = max-centering]
    tmpden = zeros(ngrid+1,1);
    for ii=1:ngrid+1
        tmpden(ii) = lp_psi(psigrid(ii));
    end
    tmpden = exp(tmpden-tmpden(psiidx));
    tmpden = tmpden/(sum(tmpden)*(psigrid(2)-psigrid(1)));
    den_psi = tmpden(psiidx);

    store_lpost(isim,:) = den_psi;
end
lpost(3) = log(mean(store_lpost));
ML = llike + lpri - sum(lpost);

out = struct('llike', llike, 'lpri', lpri, 'lpost', lpost, ...
    'store_lpost', store_lpost, 'store_lpost1', store_lpost1, ...
    'A_mean', A_mean, 'Sig_mean', Sig_mean, 'theta_mean', theta_mean);
end
