% bvar.ml.kron_bvar_csv_t_ma - log marginal likelihood of the BVAR-CSV-t-MA
% model by Chib's method: integrated likelihood by importance sampling,
% Rao-Blackwellized ordinates for (A,Sig), sigh2 and nu, and rho / psi
% ordinates from a reduced run at the posterior means. Consumes rng: R*T
% randn in the intlike, then the reduced run's draws.
%
% Two known legacy defects, both reproduced by 'bugcompat', true:
% (1) the (A,Sig) ordinate loop (legacy lines 36-51) never refreshes psi or
%     rebuilds Hpsi, so all nsims terms condition on the leftover final
%     draw; ml_BVAR_CSV_MA.m lines 26-30 show the intended per-draw pattern.
% (2) the reduced run's psi target (line 108) passes the leftover final Sig
%     where Sig_mean is intended - as ml_BVAR_CSV_MA.m line 83 does, and as
%     the h step of this same reduced run already does via CSig.
% The default path fixes both. Quirks kept verbatim in either mode (rho
% bound .999, psigrid +/-.999, ngrid 299, the lam step's missing (1+psi^2)
% correction, the reduced run's warm start) are sampler details rather than
% evaluation-point inconsistencies. tests/variant_map.md has the audit, the
% full quirk list and the effect on the published values.
%
% Extracted 2026-09-02 from
% chan2020_jbes_kronecker/legacy/ml_BVAR_CSV_t_MA.m.
%
%   [ML, out] = bvar.ml.kron_bvar_csv_t_ma(shortY, X, pri, est, ...)
%
%   pri: A0, VA0, nu0, S0, psi0, Vpsi, rho0, Vrho, nuh0, Sh0, nuub
%   est: nsims, store_A, store_Sig (running sums), store_h, store_lam,
%        store_theta ([psi rho sigh2 nu] columns), state.psihat,
%        state.psi + state.Sig (final chain draws; required under bugcompat)
%   options (name-value): 'bugcompat' (default false); 'R' - importance-
%        sampling draws (default 10000 = legacy lines 20-21); 'nsims2' -
%        reduced-run length (default 1000 = legacy line 71)
%   out: llike, lpri, lpost, store_lpost ([den_rho den_psi] reduced-run
%        columns), store_lpost1, A_mean, Sig_mean, theta_mean, bugcompat
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics, 38(1), 68-79.

function [ML, out] = kron_bvar_csv_t_ma(shortY, X, pri, est, varargin)
bugcompat = false;
R = 10000;                              % legacy ml_BVAR_CSV_t_MA.m lines 20-21
nsims2 = 1000;                          % legacy line 71
for iv = 1:2:numel(varargin)
    switch lower(varargin{iv})
        case 'bugcompat', bugcompat = varargin{iv+1};
        case 'r', R = varargin{iv+1};
        case 'nsims2', nsims2 = varargin{iv+1};
        otherwise, error('bvar:ml:kron_bvar_csv_t_ma:badOption', ...
                'unknown option ''%s''', varargin{iv});
    end
end
[T, n] = size(shortY);
k = size(X, 2);
A0 = pri.A0; VA0 = pri.VA0; nu0 = pri.nu0; S0 = pri.S0;
psi0 = pri.psi0; Vpsi = pri.Vpsi;
rho0 = pri.rho0; Vrho = pri.Vrho; nuh0 = pri.nuh0; Sh0 = pri.Sh0; nuub = pri.nuub;
nsims = est.nsims;
store_A = est.store_A; store_Sig = est.store_Sig;
store_h = est.store_h; store_lam = est.store_lam; store_theta = est.store_theta;
psihat = est.state.psihat;              % estimation run's final psi-MH mode
options = optimset('Display', 'off', 'LargeScale','off') ;  % = BVAR_CSV_t_MA.m line 31
lpri_psi = @(x) -.5*(x-psi0)^2/Vpsi -1e10*(x<-.99 || x>.99);    % = BVAR_CSV_t_MA.m line 10

ngrid = 299;                            % [ml_BVAR_CSV_t_MA.m line 10]

    % [lines 12-18]
A_mean = store_A/nsims;
Sig_mean = store_Sig/nsims;
theta_mean = mean(store_theta)';
psi_mean = theta_mean(1);
rho_mean = theta_mean(2);
sigh2_mean = theta_mean(3);
nu_mean = theta_mean(4);

    % [lines 20-28]
llike = bvar.ml.intlike_csv_t_ma(shortY,X,A_mean,Sig_mean,psi_mean,rho_mean,...
    sigh2_mean,nu_mean,R);
c_psi = 1/(normcdf(1,psi0,sqrt(Vpsi))-normcdf(-1,psi0,sqrt(Vpsi)));
c_rho = 1/(normcdf(1,rho0,sqrt(Vrho))-normcdf(-1,rho0,sqrt(Vrho)));
lpri = -.5*log(2*pi*Vpsi) + log(c_psi) -.5*(theta_mean(1)-psi0)^2/Vpsi ...
    + log(1/(nuub-2)) ...
    -.5*log(2*pi*Vrho) + log(c_rho) -.5*(rho_mean-rho0)^2/Vrho ...
    + nuh0*log(Sh0) - gammaln(nuh0) - (nuh0+1)*log(sigh2_mean) - Sh0/sigh2_mean ...
    + bvar.ml.lniwpdf(A_mean,Sig_mean,A0,sparse(1:k,1:k,1./VA0),nu0,S0);

    % evaluate the posterior density [lines 31-69]
store_lpost = zeros(nsims,3); % [log density of A Sig, log density of sigh2, log density of nu]
nugrid = sort([nu_mean; linspace(2,nuub,700)']);
nuidx = find(nugrid==nu_mean);

if bugcompat
        % Legacy defect (1) reproduced: psi and Hpsi frozen at the
        % estimation run's final chain values for the whole ordinate loop
    psi = est.state.psi;
    Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
end
for isim = 1:nsims
    h = store_h(isim,:)';
    lam = store_lam(isim,:)';
    rho = store_theta(isim,2);
    if ~bugcompat
            % corrected: integrate over the stored psi draws (the model-7
            % pattern, ml_BVAR_CSV_MA.m lines 26-30)
        psi = store_theta(isim,1);
        Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
    end

        % compute the conditional density of Sig and A
    Xtld = Hpsi\X;
    Ytld = Hpsi\shortY;
    iO_h_lam_psi = sparse(1:T,1:T,[1/(1+psi^2)*exp(-h(1));exp(-h(2:end))]./lam);
    XiO = Xtld'*iO_h_lam_psi;
    KA = sparse(1:k,1:k,1./VA0) + XiO*Xtld;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiO*Ytld);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + Ytld'*iO_h_lam_psi*Ytld ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    lden_ASig = bvar.ml.lniwpdf(A_mean,Sig_mean,Ahat,KA,nu0+T,Shat);

        % compute the conditional density of sigh2
    eh = [h(1)*sqrt(1-rho^2);  h(2:end)-rho*h(1:end-1)];
    lden_sigh2 = bvar.ml.linvgammpdf(sigh2_mean,nuh0+T/2,Sh0 + sum(eh.^2)/2);

        % compute the conditional density of nu
    sum1 = sum(log(lam));
    sum2 = sum(1./lam);
    fnu = @(x) T*(x/2.*log(x/2)-gammaln(x/2)) - (x/2+1)*sum1 - x/2*sum2;
    tmpden = fnu(nugrid);
    tmpden = exp(tmpden-max(tmpden));
    tmpden = tmpden/(sum(tmpden)*(nugrid(2)-nugrid(1)));
    den_nu = log(tmpden(nuidx));

    store_lpost(isim,:) = [lden_ASig lden_sigh2 den_nu];
end
tmpmax = max(store_lpost);
lpost = log(mean(exp(store_lpost-repmat(tmpmax,nsims,1)))) + tmpmax;

    % rho and psi ordinates: reduced run [lines 71-149]; (h, lam, rho)
    % continue from the last stored draws, psi from the final chain value
    % (= store_theta(nsims,1) in both modes), psihat from the estimation run
store_lpost1 = store_lpost;
store_lpost = zeros(nsims2,2); % [log density of rho, log density of rho of psi]
rhogrid = sort([rho_mean; linspace(-.999,.999,ngrid)']);
rhoidx = find(rhogrid==rho_mean);
psigrid = sort([theta_mean(1); linspace(-.999,.999,ngrid)']);
psiidx = find(psigrid==theta_mean(1));
U = shortY - X*A_mean;
CSig = chol(Sig_mean,'lower');
nu = nu_mean;
sigh2 = sigh2_mean;
if bugcompat
        % Legacy defect (2) reproduced: the psi target conditions on the
        % final estimation DRAW of Sig, not Sig_mean
    Sig_psi = est.state.Sig;
else
    Sig_psi = Sig_mean;
end
if ~bugcompat
        % entry state after the corrected ordinate loop: psi/Hpsi carry the
        % last stored draw - the same values the bugcompat path froze
    psi = store_theta(nsims,1);
    Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
end

for isim = 1:nsims2
    Utld = Hpsi\U;
    tmp = (Utld/CSig');
    s2_h = sum(tmp.^2,2)./lam;
    s2_h(1) = s2_h(1)/(1+psi^2);
    h = bvar.sv.csv_armh(s2_h,rho,sigh2,h,n);            % legacy line 87: root sample_h

        % sample lam [no (1+psi^2) first-observation correction - verbatim
        % legacy line 90, mirroring the estimation script]
    s2_lam = sum(tmp.^2,2)./exp(h);
    lam = 1./gamrnd((n+nu)/2,2./(s2_lam+nu));

        % sample rho [reduced-run bound .999; model-8 estimation uses .99]
    Krho = 1/Vrho + sum(h(1:T-1).^2)/sigh2;
    rhohat = Krho\(rho0/Vrho + h(1:T-1)'*h(2:T)/sigh2);
    rhoc = rhohat + sqrt(Krho)'\randn;
    grho = @(x) -.5*log(sigh2./(1-x.^2))-.5*(1-x.^2)/sigh2*h(1)^2;
    if abs(rhoc)<.999
        alpMH = exp(grho(rhoc)-grho(rho));
        if alpMH>rand
            rho = rhoc;
        end
    end

        % sample psi
    U_psi = U./repmat(sqrt(lam),1,n);
    lp_psi = @(x) bvar.ml.llike_csv_ma(x,U_psi,Sig_psi,h) + lpri_psi(x);
    if (mod(isim,100)==0) || isim == 1 %% get the Hessian every 100 iterations
        [psihat,fval,exitflag,output,grad,hess] ...
            = fminunc(@(x)-lp_psi(x),psihat,options);   %#ok<ASGLU>
        [tmpCpsi,flag] = chol(hess,'lower');            %#ok<ASGLU>
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
        Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
    end

        % compute the conditional density of rho
    tmpden = grho(rhogrid) + -.5*Krho*(rhogrid-rhohat).^2;
    tmpden = exp(tmpden-max(tmpden));
    tmpden = tmpden/(sum(tmpden)*(rhogrid(2)-rhogrid(1)));
    den_rho = tmpden(rhoidx);

        % compute the conditional density of psi
    tmpden = zeros(ngrid+1,1);
    for ii=1:ngrid+1
        tmpden(ii) = lp_psi(psigrid(ii));
    end
    tmpden = exp(tmpden-max(tmpden));
    tmpden = tmpden/(sum(tmpden)*(psigrid(2)-psigrid(1)));
    den_psi = tmpden(psiidx);

    store_lpost(isim,:) = [den_rho den_psi];
end
lpost(4:5) = log(mean(store_lpost));

ML = llike + lpri - sum(lpost);

out = struct('llike', llike, 'lpri', lpri, 'lpost', lpost, ...
    'store_lpost', store_lpost, 'store_lpost1', store_lpost1, ...
    'A_mean', A_mean, 'Sig_mean', Sig_mean, 'theta_mean', theta_mean, ...
    'bugcompat', bugcompat);
end
