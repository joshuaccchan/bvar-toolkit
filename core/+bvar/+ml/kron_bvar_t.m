% bvar.ml.kron_bvar_t - log marginal likelihood of the BVAR-t model of
% Chan (2020, JBES), Chib-style: analytic Student-t likelihood at the
% posterior means (A_mean, Sig_mean, nu_mean), NIW prior ordinate plus the
% flat nu prior 1/(nuub-2), and posterior ordinates Rao-Blackwellized over
% the stored lam draws - the (A,Sig) ordinate as a log-mean-exp of
% conditional NIW densities, the nu ordinate as the mean of grid-normalized
% conditional densities. Deterministic given the stores: consumes NO rng.
%
% Extracted 2026-09-02 (step 8, Kronecker family pass). Canonical source:
%   chan2020_jbes_kronecker/legacy/ml_BVAR_t.m  (body verbatim; lniwpdf ->
%   bvar.ml.lniwpdf). CLEAN BILL (step-8 audit): every ordinate is evaluated
% at the same (A_mean, Sig_mean, nu_mean); no leftover-workspace reads beyond
% the stores and the priors, so no bugcompat flag is needed. (Inherent legacy
% quirk kept verbatim: the nu grid normalization divides by
% nugrid(2)-nugrid(1) although inserting nu_mean makes the sorted grid
% non-uniform in one interval - shared by every grid ordinate in the family.)
%
%   [ML, out] = bvar.ml.kron_bvar_t(shortY, X, pri, est)
%
%   pri: A0, VA0, nu0, S0, nuub          [replication preset, cited there]
%   est: nsims, store_A, store_Sig (running sums), store_nu, store_lam
%   out: llike, lpri, lpost, store_lpost, A_mean, Sig_mean, nu_mean
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics, 38(1), 68-79.

function [ML, out] = kron_bvar_t(shortY, X, pri, est)
[T, n] = size(shortY);
k = size(X, 2);
A0 = pri.A0; VA0 = pri.VA0; nu0 = pri.nu0; S0 = pri.S0; nuub = pri.nuub;
nsims = est.nsims;
store_A = est.store_A; store_Sig = est.store_Sig;
store_nu = est.store_nu; store_lam = est.store_lam;

    % [ml_BVAR_t.m lines 10-12]
A_mean = store_A/nsims;
Sig_mean = store_Sig/nsims;
nu_mean = mean(store_nu)';

    % evaluate the log likelihood [lines 15-20]
CSig = chol(Sig_mean,'lower');
tmp = (shortY-X*A_mean)/CSig';
s2 = sum(tmp.^2,2);
llike = T*(gammaln((nu_mean+n)/2) - gammaln(nu_mean/2) - n/2*log(nu_mean*pi))...
    - T*sum(log(diag(CSig))) -(nu_mean+n)/2*sum(log(1+s2/nu_mean));
lpri = bvar.ml.lniwpdf(A_mean,Sig_mean,A0,sparse(1:k,1:k,1./VA0),nu0,S0) + log(1/(nuub-2));

    % evaluate the posterior density [lines 23-55]
store_lpost = zeros(nsims,2); % [log density of A Sig, density of nu]
nugrid = sort([nu_mean; linspace(2,nuub,700)']);
nuidx = find(nugrid==nu_mean);

for isim = 1:nsims
    lam = store_lam(isim,:)';

        % compute the conditional density of Sig and A
    iOm = sparse(1:T,1:T,1./lam);
    XiOm = X'*iOm;
    KA = sparse(1:k,1:k,1./VA0) + XiOm*X;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiOm*shortY);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + shortY'*iOm*shortY ...
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
lpost(1) = log(mean(exp(store_lpost(:,1)-tmpmax))) + tmpmax;
lpost(2) = log(mean(store_lpost(:,2)));
ML = llike + lpri - sum(lpost);

out = struct('llike', llike, 'lpri', lpri, 'lpost', lpost, ...
    'store_lpost', store_lpost, 'A_mean', A_mean, 'Sig_mean', Sig_mean, ...
    'nu_mean', nu_mean);
end
