% bvar.ml.kron_bvar_t_csv - log marginal likelihood of the BVAR-t-CSV model of
% Chan (2020, JBES): integrated likelihood by importance sampling
% (bvar.ml.intlike_t_csv, lam analytic and h by IS) at the posterior means,
% prior ordinate, posterior ordinates for (A,Sig), sigh2 and nu
% Rao-Blackwellized over the stored (h, lam, rho) draws, and a rho ordinate
% from a reduced run of nsims sweeps re-drawing (h, lam, rho) at fixed
% (A_mean, Sig_mean, nu_mean, sigh2_mean), continuing from the final stored
% draws. Consumes rng: R*T randn in the intlike, then the reduced run's
% draws.
%
% Extracted 2026-09-02 from chan2020_jbes_kronecker/legacy/ml_BVAR_t_CSV.m,
% body verbatim. Clean bill: every ordinate sits at the same starred point;
% the leftover-workspace reads are chain continuation from the last stored
% (h, lam, rho), made explicit here. No bugcompat flag needed. Verbatim
% quirks (reduced-run rho bound .999 vs the estimation script's .9999) are
% listed in tests/variant_map.md.
%
%   [ML, out] = bvar.ml.kron_bvar_t_csv(shortY, X, pri, est, ...)
%
%   pri: A0, VA0, nu0, S0, rho0, Vrho, nuh0, Sh0, nuub
%   est: nsims, store_A, store_Sig (running sums), store_h, store_lam,
%        store_theta ([nu rho sigh2] columns)
%   options (name-value): 'R' - importance-sampling draws (default 1000 =
%        legacy ml_BVAR_t_CSV.m lines 13-14)
%   out: llike, lpri, lpost, store_lpost (reduced-run den_rho column),
%        store_lpost1, A_mean, Sig_mean, theta_mean, h_mean
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics, 38(1), 68-79.

function [ML, out] = kron_bvar_t_csv(shortY, X, pri, est, varargin)
R = 1000;                               % legacy ml_BVAR_t_CSV.m lines 13-14
for iv = 1:2:numel(varargin)
    switch lower(varargin{iv})
        case 'r', R = varargin{iv+1};
        otherwise, error('bvar:ml:kron_bvar_t_csv:badOption', ...
                'unknown option ''%s''', varargin{iv});
    end
end
[T, n] = size(shortY);
k = size(X, 2);
A0 = pri.A0; VA0 = pri.VA0; nu0 = pri.nu0; S0 = pri.S0;
rho0 = pri.rho0; Vrho = pri.Vrho; nuh0 = pri.nuh0; Sh0 = pri.Sh0; nuub = pri.nuub;
nsims = est.nsims;
store_A = est.store_A; store_Sig = est.store_Sig;
store_h = est.store_h; store_lam = est.store_lam; store_theta = est.store_theta;

    % estimation-tail posterior means [BVAR_t_CSV.m lines 121-124]
A_mean = store_A/nsims;
Sig_mean = store_Sig/nsims;
theta_mean = mean(store_theta)';
h_mean = mean(store_h)';                                %#ok<NASGU> % legacy tail value, unused by this ml script

    % [ml_BVAR_t_CSV.m lines 10-19]
nu_mean = theta_mean(1);
rho_mean = theta_mean(2);
sigh2_mean = theta_mean(3);
llike = bvar.ml.intlike_t_csv(shortY,X,A_mean,Sig_mean,rho_mean,sigh2_mean,...
    nu_mean,R);
c_rho = 1/(normcdf(1,rho0,sqrt(Vrho))-normcdf(-1,rho0,sqrt(Vrho)));
lpri = log(1/(nuub-2)) ...
    -.5*log(2*pi*Vrho) + log(c_rho) -.5*(rho_mean-rho0)^2/Vrho ...
    + nuh0*log(Sh0) - gammaln(nuh0) - (nuh0+1)*log(sigh2_mean) - Sh0/sigh2_mean ...
    + bvar.ml.lniwpdf(A_mean,Sig_mean,A0,sparse(1:k,1:k,1./VA0),nu0,S0);

    % evaluate the posterior density [lines 22-57]
store_lpost = zeros(nsims,3); % [log density of A Sig, log density of sigh2, log density of nu]
nugrid = sort([nu_mean; linspace(2,nuub,700)']);
nuidx = find(nugrid==nu_mean);

for isim = 1:nsims
    h = store_h(isim,:)';
    lam = store_lam(isim,:)';
    rho = store_theta(isim,2);

        % compute the conditional density of Sig and A
    iOm = sparse(1:T,1:T,exp(-h)./lam);
    XiOm = X'*iOm;
    KA = sparse(1:k,1:k,1./VA0) + XiOm*X;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiOm*shortY);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + shortY'*iOm*shortY ...
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

    % rho ordinate: reduced run [lines 59-95]; chain continuation from the
    % last stored draws (the values the legacy first loop leaves behind)
store_lpost1 = store_lpost;
store_lpost = zeros(nsims,1); % [log density of rho]
rhogrid = sort([rho_mean; linspace(-.999,.999,700)']);
rhoidx = find(rhogrid==rho_mean);
U = shortY - X*A_mean;
CSig = chol(Sig_mean,'lower');
tmp = U/CSig';
s2 = sum(tmp.^2,2);                                     %#ok<NASGU> % legacy line 65: immediately overwritten inside the loop
nu = nu_mean;
sigh2 = sigh2_mean;
for isim = 1:nsims
    s2 = sum(tmp.^2,2)./lam;
    h = bvar.sv.csv_armh(s2,rho,sigh2,h,n);              % legacy line 70: root sample_h

        % sample lam
    s2 = sum(tmp.^2,2)./exp(h);
    lam = 1./gamrnd((n+nu)/2,2./(s2+nu));

        % sample rho [reduced-run bound .999; estimation uses .9999]
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

        % compute the conditional density of rho
    tmpden = grho(rhogrid) + -.5*Krho*(rhogrid-rhohat).^2;
    tmpden = exp(tmpden-max(tmpden));
    tmpden = tmpden/(sum(tmpden)*(rhogrid(2)-rhogrid(1)));
    den_rho = tmpden(rhoidx);
    store_lpost(isim,:) = den_rho;
end
lpost(4) = log(mean(store_lpost));
ML = llike + lpri - sum(lpost);

out = struct('llike', llike, 'lpri', lpri, 'lpost', lpost, ...
    'store_lpost', store_lpost, 'store_lpost1', store_lpost1, ...
    'A_mean', A_mean, 'Sig_mean', Sig_mean, 'theta_mean', theta_mean);
end
