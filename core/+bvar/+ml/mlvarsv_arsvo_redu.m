% bvar.ml.mlvarsv_arsvo_redu - log marginal likelihood of the VAR-SVO model
% (Cholesky SV with an outlier component) by adaptive importance sampling. Same
% construction as bvar.ml.mlvarsv_arsv_redu, plus the outlier scales o_t and the
% outlier probability po: o_t is drawn from its smoothed empirical posterior
% (o_hat, one categorical distribution per period over the 32 grid atoms) and po
% from a fitted beta.
%
% Three legacy defects, all in the o block, all reproduced by 'bugcompat', true:
% (1) legacy line 180 spreads the outlier prior mass over numel(o_grid) = 32
%     atoms where the sampler spreads it over 31 (VAR_ARSVO_redu.m 8/112);
% (2) line 181 indexes the T x 32 o_hat linearly, so for T >= 32 every index
%     lands in column 1 instead of at (t, o_idx(t));
% (3) line 153 omits -n*sum(log(o)), the Jacobian of the outlier scaling that
%     line 155 applies and the sampler's o step carries (VAR_ARSVO_redu.m 115).
% The default path fixes all three. None consumes rng, so both modes draw the
% identical stream and differ only in the weights. tests/variant_map.md has the
% audit, the effect on the published value and the family-wide quirks.
%
% Extracted 2026-09-03 from chan2023_joe_mlvarsv/legacy/utility/ml_var_arsvo_redu.m
% (body verbatim apart from the three corrections; helper calls redirected to core).
%
%   [lml,lmlstd,out] = bvar.ml.mlvarsv_arsvo_redu(X,Y,Y0,M,Hyper,flag_marg,...
%       store_h,store_beta,store_hpara,store_kappa,store_o,store_po,o_grid,...
%       is_kappafixed,is_kappasym, 'bugcompat',false)
%
%   flag_marg - 2 only (as in the legacy switch)
%   Hyper: as bvar.ml.mlvarsv_arsv_redu plus p0a, p0b (the beta prior on po)
%   store_o  - nsim x T outlier scales      [VAR_ARSVO_redu.m 15]
%   store_po - nsim x 1                     [16]
%   o_grid   - the 32-atom grid [1; linspace(2,20,31)'] the sampler used [9].
%              o_hat is built by exact == against store_o, so this must be the
%              same grid the chain drew from (run_ml passes out.o_grid)
%   out: store_w, bigml, store_lr_o, store_lJ_o (the o Jacobian actually
%        applied; all zeros under bugcompat), o_hat, and the fitted IS parameters
%
% Core used: bvar.priors.minn (legacy prior_Minn, n0pre = 4), bvar.priors.impact_B0
% (prior_B0), bvar.util.tnormrnd, bvar.util.vec, bvar.util.ldet, bvar.ml.isden_arss.
%
% See:
% Chan, J.C.C. (2023). Comparing stochastic volatility specifications for large
% Bayesian VARs, Journal of Econometrics, 235(2), 1419-1446.

function [lml,lmlstd,out] = mlvarsv_arsvo_redu(X,Y,Y0,M,Hyper,flag_marg,store_h,...
    store_beta,store_hpara,store_kappa,store_o,store_po,o_grid,is_kappafixed,is_kappasym,varargin)
bugcompat = false;
for iv = 1:2:numel(varargin)
    switch lower(varargin{iv})
        case 'bugcompat', bugcompat = varargin{iv+1};
        otherwise, error('bvar:ml:mlvarsv_arsvo_redu:badOption', ...
                'unknown option ''%s''', varargin{iv});
    end
end
assert(isequal(flag_marg,2), 'bvar:ml:mlvarsv_arsvo_redu:flagMarg', ...
    'only flag_marg = 2 is implemented (as in the legacy switch)');
[T,n] = size(Y);
k = size(X,2);
p = (k-1)/n;
k_beta = n*(n-1)/2;       % dimension of B0
ngrid = size(o_grid,1);   % number of atoms (32); the sampler ngrid is one fewer
M = 50*ceil(M/50);
kappa3 = 100;             % [ml_var_arsvo_redu.m 14]
B0_id = nonzeros(tril(reshape(1:n^2,n,n),-1)');
B0 = eye(n);

    % how the prior mass po is split, and how o_hat is indexed (defects 1-2)
if bugcompat
    n_pri = ngrid;                        % legacy: po/32, o_lpri a 33-vector
else
    n_pri = ngrid - 1;                    % the sampler's 31 outlier atoms
end
t_idx = (1:T)';

   % obtain parameters for importance sampling densities
if is_kappafixed % kappa1 and kappa2 are fixed, but kappa4 is free
    tmp_hat = gamfit(store_kappa(:,3));
    big_kappa4 = gamrnd(tmp_hat(1),tmp_hat(2),M,1);
    ckappa_hat = [tmp_hat(1); 1/tmp_hat(2)];
elseif is_kappasym % kappa1 = kappa2, kappa4 is free
    ckappa_hat = zeros(2,2);
    big_kappa = zeros(M,2);
    for i=2:3
        tmp_hat = gamfit(store_kappa(:,i));
        big_kappa(:,i-1) = gamrnd(tmp_hat(1),tmp_hat(2),M,1);
        ckappa_hat(i-1,:) = [tmp_hat(1); 1/tmp_hat(2)];
    end
else
    ckappa_hat = zeros(3,2);
    big_kappa = zeros(M,3);
    for i=1:3
        tmp_hat = gamfit(store_kappa(:,i));
        big_kappa(:,i) = gamrnd(tmp_hat(1),tmp_hat(2),M,1);
        ckappa_hat(i,:) = [tmp_hat(1); 1/tmp_hat(2)];
    end
end
betahat = mean(store_beta)';
betacov = cov(store_beta);
betaCcov = sparse(chol(betacov,'lower'));
muhat = mean(store_hpara(:,1:n))';
muvar = var(store_hpara(:,1:n))';
phihat = mean(store_hpara(:,n+1:2*n))';
phivar = var(store_hpara(:,n+1:2*n))';
tmp = zeros(n,2);
for i=1:n
    tmp(i,:) = gamfit(1./store_hpara(:,2*n+i));   % dead: the lines that read it are commented out [51]
end
% nusig2hat = tmp(:,1); Ssig2hat = 1./tmp(:,2);
h_hat = zeros(T*n,1);
for ii=1:n
    [hi_hat,Khi_hat] = bvar.ml.isden_arss(store_h(:,:,ii));
    h_hat((ii-1)*T+1:ii*T) = hi_hat;
    tmp0 = sparse((ii-1)*T,T);
    if ii == 1
        Kh_hat = Khi_hat;
    else
        Kh_hat = [[Kh_hat tmp0]; [tmp0' Khi_hat]];   % block-diagonal, one series at a time
    end
end
CKh_hat = chol(Kh_hat,'lower');
tmp_count = zeros(T,ngrid);
for ii = 1:ngrid
    tmp_count(:,ii) = sum(store_o == o_grid(ii));
end
o_hat = (tmp_count+1)/(size(store_o,1) + ngrid);
cumsumo_hat = cumsum(o_hat,2);
po_hat = betafit(store_po);

    % obtain importance sampling draws
big_mu = repmat(muhat',M,1) + repmat(sqrt(muvar)',M,1).*randn(M,n);
big_phi = zeros(M,n);
for ii=1:n
   big_phi(:,ii) = bvar.util.tnormrnd(phihat(ii),phivar(ii),-1,1,M);
end
% big_sig2 = zeros(M,n);
% for i=1:n
%     big_sig2(:,i) = 1./gamrnd(nusig2hat(i),1./Ssig2hat(i),M,1);
% end
big_po = betarnd(po_hat(1),po_hat(2),M,1);

switch flag_marg
    case 2
        if is_kappafixed
            cprior = -.5*(2*n)*log(2*pi) -sum(log(normcdf((1-Hyper.phi0)./sqrt(Hyper.Vphi)) -normcdf((-1-Hyper.phi0)./sqrt(Hyper.Vphi)))) ...
                -.5*sum(log(Hyper.Vmu)) -.5*sum(log(Hyper.Vphi))...
                +Hyper.c0(3,1)*log(Hyper.c0(3,2)) -gammaln(Hyper.c0(3,1));
            prior = @(m,ph,xk4,xpo) cprior -.5*sum((m-Hyper.mu0).^2./Hyper.Vmu)...
                -.5*sum((ph-Hyper.phi0).^2./Hyper.Vphi)...
                +(Hyper.c0(3,1)-1)*log(xk4) -Hyper.c0(3,2)'*xk4 ...
                + log(betapdf(xpo, Hyper.p0a, Hyper.p0b));
            cIS = -.5*(2*n)*log(2*pi) -sum(log(normcdf((1-phihat)./sqrt(phivar)) -normcdf((-1-phihat)./sqrt(phivar))))...
                -.5*sum(log(muvar)) -.5*sum(log(phivar))...
                +ckappa_hat(1)*log(ckappa_hat(2)) -gammaln(ckappa_hat(1));
            gIS = @(m,ph,xk4,xpo) cIS -.5*sum((m-muhat).^2./muvar) -.5*sum((ph-phihat).^2./phivar)...
                +(ckappa_hat(1)-1)'*log(xk4) -ckappa_hat(2)'*xk4 ...
                + log(betapdf(xpo, po_hat(1), po_hat(2)));
        elseif is_kappasym
            cprior = -.5*(2*n)*log(2*pi) -sum(log(normcdf((1-Hyper.phi0)./sqrt(Hyper.Vphi)) -normcdf((-1-Hyper.phi0)./sqrt(Hyper.Vphi)))) ...
                -.5*sum(log(Hyper.Vmu)) -.5*sum(log(Hyper.Vphi)) ...
                + Hyper.c0(2:3,1)'*log(Hyper.c0(2:3,2)) - sum(gammaln(Hyper.c0(2:3,1)));
            prior = @(m,ph,xk,xpo) cprior -.5*sum((m-Hyper.mu0).^2./Hyper.Vmu) -.5*sum((ph-Hyper.phi0).^2./Hyper.Vphi) ...
                +(Hyper.c0(2:3,1)-1)'*log(xk) - Hyper.c0(2:3,2)'*xk ...
                +log(betapdf(xpo, Hyper.p0a, Hyper.p0b));
            cIS = -.5*(2*n)*log(2*pi) -sum(log(normcdf((1-phihat)./sqrt(phivar)) -normcdf((-1-phihat)./sqrt(phivar)))) ...
                -.5*sum(log(muvar)) -.5*sum(log(phivar)) ...
                +ckappa_hat(:,1)'*log(ckappa_hat(:,2)) -sum(gammaln(ckappa_hat(:,1)));
            gIS = @(m,ph,xk,xpo) cIS -.5*sum((m-muhat).^2./muvar) -.5*sum((ph-phihat).^2./phivar)...
                +(ckappa_hat(:,1)-1)'*log(xk) -ckappa_hat(:,2)'*xk ...
                +log(betapdf(xpo, po_hat(1), po_hat(2)));
        else
            cprior = -.5*(2*n)*log(2*pi) -sum(log(normcdf((1-Hyper.phi0)./sqrt(Hyper.Vphi)) -normcdf((-1-Hyper.phi0)./sqrt(Hyper.Vphi)))) ...
                -.5*sum(log(Hyper.Vmu)) -.5*sum(log(Hyper.Vphi)) ...
                +Hyper.c0(:,1)'*log(Hyper.c0(:,2)) - sum(gammaln(Hyper.c0(:,1)));
            prior = @(m,ph,xk,xpo) cprior -.5*sum((m-Hyper.mu0).^2./Hyper.Vmu) -.5*sum((ph-Hyper.phi0).^2./Hyper.Vphi) ...
                +(Hyper.c0(:,1)-1)'*log(xk) -Hyper.c0(:,2)'*xk ...
                + log(betapdf(xpo, Hyper.p0a, Hyper.p0b));
            cIS = -.5*(2*n)*log(2*pi) -sum(log(normcdf((1-phihat)./sqrt(phivar)) -normcdf((-1-phihat)./sqrt(phivar)))) ...
                -.5*sum(log(muvar)) -.5*sum(log(phivar)) ...
                +ckappa_hat(:,1)'*log(ckappa_hat(:,2)) -sum(gammaln(ckappa_hat(:,1)));
            gIS = @(m,ph,xk,xpo) cIS -.5*sum((m-muhat).^2./muvar) -.5*sum((ph-phihat).^2./phivar)...
                +(ckappa_hat(:,1)-1)'*log(xk) -ckappa_hat(:,2)'*xk...
                +log(betapdf(xpo, po_hat(1), po_hat(2)));
        end
end
store_w = zeros(M,1);
store_lr_o = zeros(M,1);
store_lJ_o = zeros(M,1);
c_h = -T*n/2*log(2*pi) + .5*bvar.util.ldet(Kh_hat);
for isim = 1:M
    if is_kappafixed
       kappa = [store_kappa(1,1:2),big_kappa4(isim)]';
    elseif is_kappasym
       kappa = [big_kappa(isim,1),big_kappa(isim,1),big_kappa(isim,2)]';
    else
       kappa = big_kappa(isim,:)';
    end
    kappa1 = kappa(1);
    kappa2 = kappa(2);
    kappa4 = kappa(3);
    [Hyper.alp0,Hyper.Valp] = bvar.priors.minn(p,kappa1,kappa2,kappa3,Y0,Y,4);
    [Hyper.beta0,Hyper.Vbeta] = bvar.priors.impact_B0(Y0,Y,kappa4);
    mu = big_mu(isim,:)';
    phi = big_phi(isim,:)';
    longh = h_hat + CKh_hat'\randn(T*n,1);
    h = reshape(longh,T,n);
    beta = betahat + betaCcov*randn(k_beta,1);
    B0(B0_id) = beta;
    po = big_po(isim);
    o_idx = ngrid - sum(repmat(rand(T,1),1,ngrid)<cumsumo_hat,2)+1;
    o = o_grid(o_idx);
        % defect 3: the legacy leaves the outlier Jacobian out of c1
    if bugcompat
        lJ_o = 0;
    else
        lJ_o = -n*sum(log(o));
    end

    c1 = -n*T/2*log(2*pi) -.5*sum(sum(h)) -.5*sum(log(Hyper.Valp)) + lJ_o;
    iValp = sparse(1:n*k,1:n*k,1./Hyper.Valp);
    diag_sqrt_D = bvar.util.vec(exp(h/2).*repmat(o,1,n));
    ytilde = bvar.util.vec(Y*B0')./diag_sqrt_D;
    Xtilde = kron(B0,X)./diag_sqrt_D;
    Kalp = iValp + Xtilde'*Xtilde;
    CKalp = chol(Kalp,'lower');
    tmpc = CKalp\(iValp*Hyper.alp0 + Xtilde'*ytilde);
    lclike = c1 -sum(log(diag(CKalp))) -.5*(sum(ytilde.^2) +sum(Hyper.alp0.^2./Hyper.Valp) -tmpc'*tmpc);
    lh_pri = 0;
    switch flag_marg
        case 2
            for ii=1:n
                mui = mu(ii);
                phii = phi(ii);
                hi = h(:,ii);
                Hphii = speye(T) -phii*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
                HiSH = Hphii'*sparse(1:T,1:T,[1-phii^2, ones(1,T-1)])*Hphii;
                lh_pri = lh_pri -T/2*log(2*pi) +.5*log(1-phii^2)...
                    +Hyper.nuh(ii)*log(Hyper.Sh(ii)) -gammaln(Hyper.nuh(ii))...
                    +gammaln(Hyper.nuh(ii)+T/2) -(Hyper.nuh(ii)+T/2)*log(Hyper.Sh(ii)+.5*(hi-mui)'*HiSH*(hi-mui));
            end
    end
    lh_g = c_h -.5*(longh-h_hat)'*Kh_hat*(longh-h_hat);
    tmpbeta = betaCcov\(beta-betahat);
    lr_beta = -.5*sum(log(Hyper.Vbeta)) -.5*sum((beta-Hyper.beta0).^2./Hyper.Vbeta) ...
        +sum(log(diag(betaCcov))) +.5*(tmpbeta'*tmpbeta);
        % defects 1-2: prior mass over n_pri atoms, and how o_hat is indexed
    o_lpri = log([1-po; repmat(po/n_pri,n_pri,1)]);
    if bugcompat
        lin_o = o_idx;                    % linear index: column 1 of o_hat for T >= 32
    else
        lin_o = (o_idx-1)*T + t_idx;      % (t, o_idx(t))
    end
    lr_o = sum(o_lpri(o_idx) - log(o_hat(lin_o)));
    if is_kappafixed
        store_w(isim) = lclike + lh_pri - lh_g + lr_beta + lr_o + prior(mu,phi,kappa4,po) - gIS(mu,phi,kappa4,po);
    elseif is_kappasym
        store_w(isim) = lclike + lh_pri - lh_g + lr_beta + lr_o + prior(mu,phi,kappa(2:3),po) - gIS(mu,phi,kappa(2:3),po);
    else
        store_w(isim) = lclike + lh_pri - lh_g + lr_beta + lr_o + prior(mu,phi,kappa,po) - gIS(mu,phi,kappa,po);
    end
    store_lr_o(isim) = lr_o;
    store_lJ_o(isim) = lJ_o;
end
shortw = reshape(store_w,M/50,50);
maxw = max(shortw);
bigml = log(mean(exp(shortw-repmat(maxw,M/50,1)),1)) + maxw;
lml = mean(bigml);
lmlstd = std(bigml)/sqrt(50);

out = struct('store_w',store_w, 'bigml',bigml, 'M',M, 'bugcompat',bugcompat, ...
    'store_lr_o',store_lr_o, 'store_lJ_o',store_lJ_o, ...
    'o_hat',o_hat, 'po_hat',po_hat, 'ngrid_atoms',ngrid, 'n_prior_atoms',n_pri, ...
    'h_hat',h_hat, 'Kh_hat',Kh_hat, 'betahat',betahat, 'betacov',betacov, ...
    'muhat',muhat, 'muvar',muvar, 'phihat',phihat, 'phivar',phivar, ...
    'ckappa_hat',ckappa_hat);
end
