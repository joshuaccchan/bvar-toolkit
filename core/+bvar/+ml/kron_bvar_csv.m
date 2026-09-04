% bvar.ml.kron_bvar_csv - log marginal likelihood of the BVAR-CSV model of
% Chan (2020, JBES): integrated likelihood by importance sampling
% (bvar.ml.intlike_csv, h integrated out) at the posterior means, prior
% ordinate (NIW x truncated-normal rho x inverse-gamma sigh2), posterior
% ordinates for (A,Sig) and sigh2 Rao-Blackwellized over the stored (h, rho)
% draws, and a rho ordinate from a REDUCED MCMC RUN of nsims sweeps that
% re-draws (h, rho) at fixed (A_mean, Sig_mean, sigh2_mean), continuing the
% chain from the final stored draws. Consumes rng: R*T randn in the intlike,
% then the reduced run's AR-MH h draws and rho MH draws.
%
% Extracted 2026-09-02 from chan2020_jbes_kronecker/legacy/ml_BVAR_CSV.m,
% body verbatim. Clean bill: every ordinate sits at the same starred point;
% the leftover-workspace reads are chain continuation, made explicit here
% (reduced run starts from the last stored h and rho; est.state.countrho
% carries the estimation counter the legacy script increments). No bugcompat
% flag needed. The reduced run's h step is bvar.sv.csv_armh with NR start
% h_mean and a first-sweep forced accept.
%
%   [ML, out] = bvar.ml.kron_bvar_csv(shortY, X, pri, est, ...)
%
%   pri: A0, VA0, nu0, S0, rho0, Vrho, nuh0, Sh0
%   est: nsims, store_A, store_Sig (running sums), store_h, store_theta
%        ([rho sigh2] columns), state.countrho
%   options (name-value): 'R' - importance-sampling draws for the integrated
%        likelihood (default 1000 = legacy ml_BVAR_CSV.m line 10)
%   out: llike, lpri, lpost, store_lpost (reduced-run den_rho column),
%        countrho, A_mean, Sig_mean, theta_mean, h_mean
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics, 38(1), 68-79.

function [ML, out] = kron_bvar_csv(shortY, X, pri, est, varargin)
R = 1000;                               % legacy ml_BVAR_CSV.m line 10
for iv = 1:2:numel(varargin)
    switch lower(varargin{iv})
        case 'r', R = varargin{iv+1};
        otherwise, error('bvar:ml:kron_bvar_csv:badOption', ...
                'unknown option ''%s''', varargin{iv});
    end
end
[T, n] = size(shortY);
k = size(X, 2);
A0 = pri.A0; VA0 = pri.VA0; nu0 = pri.nu0; S0 = pri.S0;
rho0 = pri.rho0; Vrho = pri.Vrho; nuh0 = pri.nuh0; Sh0 = pri.Sh0;
nsims = est.nsims;
store_A = est.store_A; store_Sig = est.store_Sig;
store_h = est.store_h; store_theta = est.store_theta;
countrho = est.state.countrho;          % legacy: continues the estimation counter

    % estimation-tail posterior means [BVAR_CSV.m lines 95-98]
A_mean = store_A/nsims;
Sig_mean = store_Sig/nsims;
h_mean = mean(store_h)';
theta_mean = mean(store_theta)';

    % [ml_BVAR_CSV.m lines 10-14]
llike = bvar.ml.intlike_csv(shortY,X,A_mean,Sig_mean,theta_mean(1),theta_mean(2),R);
c_rho = 1/(normcdf(1,rho0,sqrt(Vrho))-normcdf(-1,rho0,sqrt(Vrho)));
lpri = -.5*log(2*pi*Vrho) + log(c_rho) -.5*(theta_mean(1)-rho0)^2/Vrho ...
    + nuh0*log(Sh0) - gammaln(nuh0) - (nuh0+1)*log(theta_mean(2)) - Sh0/theta_mean(2) ...
    + bvar.ml.lniwpdf(A_mean,Sig_mean,A0,sparse(1:k,1:k,1./VA0),nu0,S0);

    % evaluate the posterior density [lines 17-40]
store_lpost = zeros(nsims,2); % [log density of A Sig, log density of sigh2]

for isim = 1:nsims
    h = store_h(isim,:)';
    rho = store_theta(isim,1);

        % compute the conditional density of Sig and A
    iOh = sparse(1:T,1:T,exp(-h));
    XiOh = X'*iOh;
    KA = sparse(1:k,1:k,1./VA0) + XiOh*X;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiOh*shortY);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + shortY'*iOh*shortY ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    lden_ASig = bvar.ml.lniwpdf(A_mean,Sig_mean,Ahat,KA,nu0+T,Shat);

        % compute the conditional density of sigh2
    eh = [h(1)*sqrt(1-rho^2);  h(2:end)-rho*h(1:end-1)];
    lden_sigh2 = bvar.ml.linvgammpdf(theta_mean(2),nuh0+T/2,Sh0 + sum(eh.^2)/2);

    store_lpost(isim,:) = [lden_ASig lden_sigh2];
end
tmpmax = max(store_lpost);
lpost = log(mean(exp(store_lpost-repmat(tmpmax,nsims,1)))) + tmpmax;

    % rho ordinate: reduced run [lines 42-111]; chain continuation from the
    % last stored draws (the values the legacy first loop leaves behind)
store_lpost1 = store_lpost;             % first-phase densities, for out
store_lpost = zeros(nsims,1); % [log density of rho]
rhogrid = sort([theta_mean(1); linspace(-.999,.999,700)']);
rhoidx = find(rhogrid==theta_mean(1));
U = shortY - X*A_mean;
CSig = chol(Sig_mean,'lower');
tmp = (U/CSig');
s2 = sum(tmp.^2,2);
sigh2 = theta_mean(2);
for isim = 1:nsims
        % [lines 51-88] inline AR-MH h step: NR start at h_mean, forced
        % accept on the reduced run's first sweep -> bvar.sv.csv_armh
    h = bvar.sv.csv_armh(s2,rho,sigh2,h,n,isim==1,h_mean);

        % sample rho [lines 91-102]
    Krho = 1/Vrho + sum(h(1:T-1).^2)/sigh2;
    rhohat = Krho\(rho0/Vrho + h(1:T-1)'*h(2:T)/sigh2);
    rhoc = rhohat + sqrt(Krho)'\randn;
    grho = @(x) -.5*log(sigh2./(1-x.^2))-.5*(1-x.^2)/sigh2*h(1)^2;
    if abs(rhoc)<.9999
        alpMH = exp(grho(rhoc)-grho(rho));
        if alpMH>rand
            rho = rhoc;
            countrho = countrho+1;
        end
    end

        % compute the conditional density of rho [lines 105-109]
    tmpden = grho(rhogrid) + -.5*Krho*(rhogrid-rhohat).^2;
    tmpden = exp(tmpden-max(tmpden));
    tmpden = tmpden/(sum(tmpden)*(rhogrid(2)-rhogrid(1)));
    den_rho = tmpden(rhoidx);
    store_lpost(isim,:) = den_rho;
end
lpost(3) = log(mean(store_lpost));
ML = llike + lpri - sum(lpost);

out = struct('llike', llike, 'lpri', lpri, 'lpost', lpost, ...
    'store_lpost', store_lpost, 'store_lpost1', store_lpost1, ...
    'countrho', countrho, 'A_mean', A_mean, 'Sig_mean', Sig_mean, ...
    'theta_mean', theta_mean, 'h_mean', h_mean);
end
