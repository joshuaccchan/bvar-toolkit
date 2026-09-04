% bvar.ml.mlvarsv_csv - log marginal likelihood of the VAR-CSV model by adaptive
% importance sampling (not the Chib decomposition of the Kronecker family).
% (A,Sig) and the log-volatility variance sigh2 are integrated out analytically;
% h, phi and kappa are drawn from importance densities fitted to the posterior
% draws - h from bvar.ml.isden_arss, phi truncated normal, kappa gamma by
% moment/ML fit - and the M log weights are averaged in 50 batches, which also
% gives the numerical standard error.
%
% Clean bill from the step-10 audit. tests/variant_map.md has the audit and the
% family-wide quirks.
%
% Extracted 2026-09-03 from chan2023_joe_mlvarsv/legacy/utility/ml_var_csv.m
% (body verbatim; helper calls redirected to core - see below).
%
%   [lml,lmlstd,out] = bvar.ml.mlvarsv_csv(X,Y,Y0,M,Hyper,store_h,store_hpara,...
%                                         store_kappa,is_kappafixed)
%
%   Hyper: A0, VA, nu0, S0, c0, nuh, Sh, phi0, Vphi. VA is recomputed inside
%          from each kappa draw before any read, so whichever version the
%          caller passes is irrelevant (the legacy passes the final sweep's).
%   store_h     - nsim x T log-volatility draws        [VAR_CSV.m 11]
%   store_hpara - nsim x 2, columns [phi sig2]         [VAR_CSV.m 12]
%   store_kappa - nsim x 1                             [VAR_CSV.m 13]
%   out: store_w, bigml (the 50 batch values), and the fitted IS parameters
%
% Core used: bvar.priors.niw('mlvarsv_ncp') (legacy prior_NCP), bvar.util.tnormrnd,
% bvar.util.ldet, bvar.util.mgammaln, bvar.ml.isden_arss (getISden_ARSS),
% bvar.ml.lgampdf / ltnormpdf / lmvnpdf_pcn.
%
% See:
% Chan, J.C.C. (2023). Comparing stochastic volatility specifications for large
% Bayesian VARs, Journal of Econometrics, 235(2), 1419-1446.

function [lml,lmlstd,out] = mlvarsv_csv(X,Y,Y0,M,Hyper,store_h,...
    store_hpara,store_kappa,is_kappafixed)
kappa3 = 100; % overall shrinkage parameter for the intercept [ml_var_csv.m 9]
M = 50*ceil(M/50);
[T,n] = size(Y);
k = size(Hyper.A0,1);
p = (k-1)/n;

phi_hat = mean(store_hpara(:,1));
Kphi_hat = 1/var(store_hpara(:,1));
big_phi = bvar.util.tnormrnd(phi_hat,1/Kphi_hat,-1,1,M);

if ~is_kappafixed
    tmp_hat = gamfit(store_kappa);
    ckappa_hat = [tmp_hat(1); 1/tmp_hat(2)];
    big_kappa = gamrnd(ckappa_hat(1),1./ckappa_hat(2),M,1);
end

[h_hat,Kh_hat] = bvar.ml.isden_arss(store_h);
CKh_hat = chol(Kh_hat,'lower');

store_w = zeros(M,1);
if is_kappafixed
    prior = @(xr) bvar.ml.ltnormpdf(xr,Hyper.phi0,Hyper.Vphi,-1,1);
    gIS = @(xh,xr) bvar.ml.lmvnpdf_pcn(xh,h_hat,Kh_hat) + bvar.ml.ltnormpdf(xr,phi_hat,1/Kphi_hat,-1,1);
else
    prior = @(xr,xk) bvar.ml.ltnormpdf(xr,Hyper.phi0,Hyper.Vphi,-1,1) ...
        +bvar.ml.lgampdf(xk,Hyper.c0(1),Hyper.c0(2));
    gIS = @(xh,xr,xk) bvar.ml.lmvnpdf_pcn(xh,h_hat,Kh_hat) + bvar.ml.ltnormpdf(xr,phi_hat,1/Kphi_hat,-1,1) ...
        +bvar.ml.lgampdf(xk,ckappa_hat(1),ckappa_hat(2));
end
for isim = 1:M
    phi = big_phi(isim,:);
    if is_kappafixed
        kappa = store_kappa(1);
    else
        kappa = big_kappa(isim);
    end
    [~,Hyper.VA] = bvar.priors.niw(p,[kappa kappa3],Y0,Y,'mlvarsv_ncp'); % update Hyper.VA given the new kappa

    h = h_hat + CKh_hat'\randn(T,1);
    Hphi = speye(T) - phi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
    HiSH = Hphi'*sparse(1:T,1:T,[1-phi^2, ones(1,T-1)])*Hphi;
    lclike = like_VAR_CSV(Y,X,Hyper,h) -T/2*log(2*pi) +.5*log(1-phi^2)...
        +Hyper.nuh*log(Hyper.Sh) -gammaln(Hyper.nuh) + gammaln(Hyper.nuh+T/2)...
        -(Hyper.nuh+T/2)*log(Hyper.Sh +.5*h'*HiSH*h);
    if is_kappafixed
        w = lclike + prior(phi) - gIS(h,phi);
    else
        w = lclike + prior(phi,kappa) - gIS(h,phi,kappa);
    end
    store_w(isim) = w;
end
shortw = reshape(store_w,M/50,50);
maxw = max(shortw);

bigml = log(mean(exp(shortw-repmat(maxw,M/50,1)),1)) + maxw;
lml = mean(bigml);
lmlstd = std(bigml)/sqrt(50);

out = struct('store_w',store_w, 'bigml',bigml, 'M',M, ...
    'h_hat',h_hat, 'Kh_hat',Kh_hat, 'phi_hat',phi_hat, 'Kphi_hat',Kphi_hat);
if is_kappafixed
    out.ckappa_hat = [];
else
    out.ckappa_hat = ckappa_hat;
end
end

% -------------------------------------------------------------------------
function lden = like_VAR_CSV(Y,X,Hyper,h)
% log p(Y | h, kappa) with (A,Sig) integrated out under the natural-conjugate
% prior - the matric-t ordinate. [ml_var_csv.m 69-83, verbatim]
[T,n] = size(Y);
k = size(Hyper.A0,1);
iOh = sparse(1:T,1:T,exp(-h));
XiOh = X'*iOh;
K_A = sparse(1:k,1:k,1./Hyper.VA) + XiOh*X;
A_hat = K_A\(sparse(1:k,1:k,Hyper.VA)\Hyper.A0 + XiOh*Y);
S_hat = Hyper.S0 + Hyper.A0'*sparse(1:k,1:k,1./Hyper.VA)*Hyper.A0 ...
    + Y'*iOh*Y - A_hat'*K_A*A_hat;
S_hat = (S_hat+S_hat')/2; % adjust for rounding errors

lden = -n*T/2*log(pi) -n/2*sum(h) -n/2*(sum(log(Hyper.VA)) +bvar.util.ldet(K_A))...
    +Hyper.nu0/2*bvar.util.ldet(Hyper.S0) -(Hyper.nu0+T)/2*bvar.util.ldet(S_hat) ...
    +bvar.util.mgammaln(n,(Hyper.nu0+T)/2) -bvar.util.mgammaln(n,Hyper.nu0/2);
end
