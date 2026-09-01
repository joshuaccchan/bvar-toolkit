% This script computes the marginal likelihood of VAR-SV
%    
% See:
% Chan, J.C.C. (2022). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, forthcoming.

function [lml,lmlstd] = ml_var_arsv_redu(X,Y,Y0,M,Hyper,flag_marg,store_h,...
    store_beta,store_hpara,store_kappa,is_kappafixed,is_kappasym)
[T,n] = size(Y);
k = size(X,2);
p = (k-1)/n;
k_beta = n*(n-1)/2;       % dimension of B0
M = 50*ceil(M/50);
kappa3 = 100;
B0_id = nonzeros(tril(reshape(1:n^2,n,n),-1)');
B0 = eye(n);

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
    tmp(i,:) = gamfit(1./store_hpara(:,2*n+i));    
end
% nusig2hat = tmp(:,1); Ssig2hat = 1./tmp(:,2);
h_hat = zeros(T*n,1);
for ii=1:n
    [hi_hat,Khi_hat] = getISden_ARSS(store_h(:,:,ii));
    h_hat((ii-1)*T+1:ii*T) = hi_hat;
    tmp0 = sparse((ii-1)*T,T);
    if ii == 1        
        Kh_hat = Khi_hat;
    else        
        Kh_hat = [[Kh_hat tmp0]; [tmp0' Khi_hat]];
    end   
end
CKh_hat = chol(Kh_hat,'lower');

    % obtain importance sampling draws
big_mu = repmat(muhat',M,1) + repmat(sqrt(muvar)',M,1).*randn(M,n);
big_phi = zeros(M,n);
for ii=1:n
   big_phi(:,ii) = tnormrnd(phihat(ii),phivar(ii),-1,1,M);
end
% big_sig2 = zeros(M,n);
% for i=1:n
%     big_sig2(:,i) = 1./gamrnd(nusig2hat(i),1./Ssig2hat(i),M,1);
% end
        
switch flag_marg  
    case 2
        if is_kappafixed
            cprior = -.5*(2*n)*log(2*pi) -sum(log(normcdf((1-Hyper.phi0)./sqrt(Hyper.Vphi)) -normcdf((-1-Hyper.phi0)./sqrt(Hyper.Vphi)))) ...
                -.5*sum(log(Hyper.Vmu)) -.5*sum(log(Hyper.Vphi))...
                +Hyper.c0(3,1)*log(Hyper.c0(3,2)) -gammaln(Hyper.c0(3,1));
            prior = @(m,ph,xk4) cprior -.5*sum((m-Hyper.mu0).^2./Hyper.Vmu)...
                -.5*sum((ph-Hyper.phi0).^2./Hyper.Vphi)...
                +(Hyper.c0(3,1)-1)*log(xk4) -Hyper.c0(3,2)'*xk4;
            cIS = -.5*(2*n)*log(2*pi) -sum(log(normcdf((1-phihat)./sqrt(phivar)) -normcdf((-1-phihat)./sqrt(phivar))))...
                -.5*sum(log(muvar)) -.5*sum(log(phivar))...
                +ckappa_hat(1)*log(ckappa_hat(2)) -gammaln(ckappa_hat(1));
            gIS = @(m,ph,xk4) cIS -.5*sum((m-muhat).^2./muvar) -.5*sum((ph-phihat).^2./phivar)...
                +(ckappa_hat(1)-1)'*log(xk4) -ckappa_hat(2)'*xk4;
        elseif is_kappasym
            cprior = -.5*(2*n)*log(2*pi) -sum(log(normcdf((1-Hyper.phi0)./sqrt(Hyper.Vphi)) -normcdf((-1-Hyper.phi0)./sqrt(Hyper.Vphi)))) ...
                -.5*sum(log(Hyper.Vmu)) -.5*sum(log(Hyper.Vphi)) ...
                + Hyper.c0(2:3,1)'*log(Hyper.c0(2:3,2)) - sum(gammaln(Hyper.c0(2:3,1)));
            prior = @(m,ph,xk) cprior -.5*sum((m-Hyper.mu0).^2./Hyper.Vmu) -.5*sum((ph-Hyper.phi0).^2./Hyper.Vphi) ...
                +(Hyper.c0(2:3,1)-1)'*log(xk) - Hyper.c0(2:3,2)'*xk;
            cIS = -.5*(2*n)*log(2*pi) -sum(log(normcdf((1-phihat)./sqrt(phivar)) -normcdf((-1-phihat)./sqrt(phivar)))) ...
                -.5*sum(log(muvar)) -.5*sum(log(phivar)) ...
                +ckappa_hat(:,1)'*log(ckappa_hat(:,2)) -sum(gammaln(ckappa_hat(:,1)));
            gIS = @(m,ph,xk) cIS -.5*sum((m-muhat).^2./muvar) -.5*sum((ph-phihat).^2./phivar)...
                +(ckappa_hat(:,1)-1)'*log(xk) -ckappa_hat(:,2)'*xk;
        else
            cprior = -.5*(2*n)*log(2*pi) -sum(log(normcdf((1-Hyper.phi0)./sqrt(Hyper.Vphi)) -normcdf((-1-Hyper.phi0)./sqrt(Hyper.Vphi)))) ...
                -.5*sum(log(Hyper.Vmu)) -.5*sum(log(Hyper.Vphi)) ...
                +Hyper.c0(:,1)'*log(Hyper.c0(:,2)) - sum(gammaln(Hyper.c0(:,1)));
            prior = @(m,ph,xk) cprior -.5*sum((m-Hyper.mu0).^2./Hyper.Vmu) -.5*sum((ph-Hyper.phi0).^2./Hyper.Vphi) ...
                +(Hyper.c0(:,1)-1)'*log(xk) -Hyper.c0(:,2)'*xk;
            cIS = -.5*(2*n)*log(2*pi) -sum(log(normcdf((1-phihat)./sqrt(phivar)) -normcdf((-1-phihat)./sqrt(phivar)))) ...
                -.5*sum(log(muvar)) -.5*sum(log(phivar)) ...
                +ckappa_hat(:,1)'*log(ckappa_hat(:,2)) -sum(gammaln(ckappa_hat(:,1)));
            gIS = @(m,ph,xk) cIS -.5*sum((m-muhat).^2./muvar) -.5*sum((ph-phihat).^2./phivar)...
                +(ckappa_hat(:,1)-1)'*log(xk) -ckappa_hat(:,2)'*xk;
        end
end
store_w = zeros(M,1);
c_h = -T*n/2*log(2*pi) + .5*ldet(Kh_hat);
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
    [Hyper.alp0,Hyper.Valp] = prior_Minn(p,kappa1,kappa2,kappa3,Y0,Y);
    [Hyper.beta0,Hyper.Vbeta] = prior_B0(Y0,Y,kappa4);    
    mu = big_mu(isim,:)';
    phi = big_phi(isim,:)';
    longh = h_hat + CKh_hat'\randn(T*n,1);    
    h = reshape(longh,T,n);
    beta = betahat + betaCcov*randn(k_beta,1);
    B0(B0_id) = beta; 
    c1 = -n*T/2*log(2*pi) -.5*sum(sum(h)) -.5*sum(log(Hyper.Valp)); 
    iValp = sparse(1:n*k,1:n*k,1./Hyper.Valp);
    diag_sqrt_D = vec(exp(h/2));
    ytilde = vec(Y*B0')./diag_sqrt_D;    
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
    if is_kappafixed        
        store_w(isim) = lclike + lh_pri - lh_g + lr_beta + prior(mu,phi,kappa4) - gIS(mu,phi,kappa4);
    elseif is_kappasym
        store_w(isim) = lclike + lh_pri - lh_g + lr_beta + prior(mu,phi,kappa(2:3)) - gIS(mu,phi,kappa(2:3));
    else
        store_w(isim) = lclike + lh_pri - lh_g + lr_beta + prior(mu,phi,kappa) - gIS(mu,phi,kappa);
    end
        
end
shortw = reshape(store_w,M/50,50);
maxw = max(shortw);
bigml = log(mean(exp(shortw-repmat(maxw,M/50,1)),1)) + maxw;
lml = mean(bigml);
lmlstd = std(bigml)/sqrt(50);
end
  