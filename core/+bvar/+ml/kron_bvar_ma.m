% bvar.ml.kron_bvar_ma - log marginal likelihood of the BVAR-MA model by
% Chib's method. Deterministic given the stores: consumes no rng.
%
%   [ML, out] = bvar.ml.kron_bvar_ma(shortY, X, pri, est, ...)
%
%   pri: A0, VA0, nu0, S0, psi0, Vpsi   (see the replication preset)
%   est: nsims, store_A, store_Sig, store_psi, state.psi
%   options: 'bugcompat' (default false)
%
% Known legacy defect, reproduced by 'bugcompat', true: ml_BVAR_MA.m line 17
% uses the leftover final-draw psi in the first observation's variance
% correction, where psi_mean is intended and used everywhere else. The
% default path uses psi_mean throughout. tests/variant_map.md has the audit
% and the effect on the published values.
%
% Extracted 2026-09-02 from chan2020_jbes_kronecker/legacy/ml_BVAR_MA.m,
% body verbatim (lniwpdf/llike_MA -> the bvar.ml equivalents; the lpri_psi
% handle reconstructed from pri.psi0/pri.Vpsi).
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics, 38(1), 68-79.

function [ML, out] = kron_bvar_ma(shortY, X, pri, est, varargin)
bugcompat = false;
for iv = 1:2:numel(varargin)
    switch lower(varargin{iv})
        case 'bugcompat', bugcompat = varargin{iv+1};
        otherwise, error('bvar:ml:kron_bvar_ma:badOption', ...
                'unknown option ''%s''', varargin{iv});
    end
end
[T, n] = size(shortY);
k = size(X, 2);
A0 = pri.A0; VA0 = pri.VA0; nu0 = pri.nu0; S0 = pri.S0;
psi0 = pri.psi0; Vpsi = pri.Vpsi;
nsims = est.nsims;
store_A = est.store_A; store_Sig = est.store_Sig; store_psi = est.store_psi;

    % estimation-tail posterior means [BVAR_MA.m lines 112-114]
A_mean = store_A/nsims;
Sig_mean = store_Sig/nsims;
psi_mean = mean(store_psi)';

    % lpri_psi reconstructed verbatim [BVAR_MA.m line 10]
lpri_psi = @(x) -.5*(x-psi0)^2/Vpsi -1e10*(x<-.99 || x>.99);

    % evaluate the log likelihood [ml_BVAR_MA.m lines 11-17]
Hpsi = speye(T) + psi_mean*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
CSig = chol(Sig_mean,'lower');
Utld = Hpsi\(shortY-X*A_mean);
tmp = (Utld/CSig');
s2 = sum(tmp.^2,2);
if bugcompat
        % Legacy defect reproduced: line 17's `psi` is the estimation run's
        % final chain draw, not psi_mean
    psi_llike = est.state.psi;
else
    psi_llike = psi_mean;   % corrected: consistent evaluation point
end
llike =  -T*n/2*log(2*pi) - T*sum(log(diag(CSig))) -n/2*log(1+psi_mean^2) ...
    -.5*(s2(1)/(1+psi_llike^2) + sum(s2(2:end)));
c_psi = 1/(normcdf(1,psi0,sqrt(Vpsi))-normcdf(-1,psi0,sqrt(Vpsi)));
lpri = bvar.ml.lniwpdf(A_mean,Sig_mean,A0,sparse(1:k,1:k,1./VA0),nu0,S0) ...
    -.5*log(2*pi*Vpsi) + log(c_psi) -.5*(psi_mean(1)-psi0)^2/Vpsi;

    % evaluate the posterior density [lines 23-44]
store_lpost = zeros(nsims,1); % [log density of A Sig]

for isim = 1:nsims
    psi = store_psi(isim,:)';
    Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);

        % compute the conditional density of Sig and A
    Xtld = Hpsi\X;
    Ytld = Hpsi\shortY;
    iO = sparse(1:T,1:T,[1/(1+psi^2) ones(1,T-1)]);
    XtldiO = Xtld'*iO;
    KA = sparse(1:k,1:k,1./VA0) + XtldiO*Xtld;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XtldiO*Ytld);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + Ytld'*iO*Ytld ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    lden_ASig = bvar.ml.lniwpdf(A_mean,Sig_mean,Ahat,KA,nu0+T,Shat);

    store_lpost(isim,:) = lden_ASig;
end
tmpmax = max(store_lpost(:,1));
lpost = log(mean(exp(store_lpost(:,1)-tmpmax))) + tmpmax;

    % psi ordinate: grid-normalized conditional density [lines 46-57]
U = shortY - X*A_mean;
lp_psi = @(x) bvar.ml.llike_ma(x,U,Sig_mean) + lpri_psi(x);
psigrid = sort([psi_mean; linspace(-.99,.99,700)']);
psiidx = find(psigrid==psi_mean);
tmpden = zeros(701,1);
for ii=1:701
    tmpden(ii) = lp_psi(psigrid(ii));
end
tmpden = exp(tmpden-max(tmpden));
tmpden = tmpden/(sum(tmpden)*(psigrid(2)-psigrid(1)));
den_psi = tmpden(psiidx);
lpost(2) = log(den_psi);
ML = llike + lpri - sum(lpost);

out = struct('llike', llike, 'lpri', lpri, 'lpost', lpost, ...
    'store_lpost', store_lpost, 'den_psi', den_psi, 'A_mean', A_mean, ...
    'Sig_mean', Sig_mean, 'psi_mean', psi_mean, 'psi_llike', psi_llike, ...
    'bugcompat', bugcompat);
end
