% This function computes the marginal likelihood of VAR-CSV
%    
% See:
% Chan, J.C.C. (2022). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, forthcoming.

function [ml,mlstd,store_w] = ml_var_csv(X,Y,Y0,M,Hyper,store_h,...
    store_hpara,store_kappa,is_kappafixed)
kappa3 = 100; % overall shrinkage parameter for the intercept
M = 50*ceil(M/50);
[T,n] = size(Y);
k = size(Hyper.A0,1);
p = (k-1)/n;
    
phi_hat = mean(store_hpara(:,1));
Kphi_hat = 1/var(store_hpara(:,1));
big_phi = tnormrnd(phi_hat,1/Kphi_hat,-1,1,M);

if ~is_kappafixed
    tmp_hat = gamfit(store_kappa);
    ckappa_hat = [tmp_hat(1); 1/tmp_hat(2)];
    big_kappa = gamrnd(ckappa_hat(1),1./ckappa_hat(2),M,1);
end

[h_hat,Kh_hat] = getISden_ARSS(store_h);
CKh_hat = chol(Kh_hat,'lower');

store_w = zeros(M,1);
if is_kappafixed
    prior = @(xr) ltnormpdf(xr,Hyper.phi0,Hyper.Vphi,-1,1);
    gIS = @(xh,xr) lmvnpdf_pcn(xh,h_hat,Kh_hat) + ltnormpdf(xr,phi_hat,1/Kphi_hat,-1,1);
else
    prior = @(xr,xk) ltnormpdf(xr,Hyper.phi0,Hyper.Vphi,-1,1) ...
        +lgampdf(xk,Hyper.c0(1),Hyper.c0(2));
    gIS = @(xh,xr,xk) lmvnpdf_pcn(xh,h_hat,Kh_hat) + ltnormpdf(xr,phi_hat,1/Kphi_hat,-1,1) ...
        +lgampdf(xk,ckappa_hat(1),ckappa_hat(2));    
end
for isim = 1:M
    phi = big_phi(isim,:); 
    if is_kappafixed
        kappa = store_kappa(1);
    else
        kappa = big_kappa(isim);
    end
    [~,Hyper.VA] = prior_NCP(p,kappa,kappa3,Y0,Y); % update Hyper.VA given the new kappa
    
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
ml = mean(bigml);
mlstd = std(bigml)/sqrt(50);

end

function lden = like_VAR_CSV(Y,X,Hyper,h)
[T,n] = size(Y);
k = size(Hyper.A0,1);
iOh = sparse(1:T,1:T,exp(-h));
XiOh = X'*iOh;
K_A = sparse(1:k,1:k,1./Hyper.VA) + XiOh*X;
A_hat = K_A\(sparse(1:k,1:k,Hyper.VA)\Hyper.A0 + XiOh*Y);
S_hat = Hyper.S0 + Hyper.A0'*sparse(1:k,1:k,1./Hyper.VA)*Hyper.A0 ...
    + Y'*iOh*Y - A_hat'*K_A*A_hat;
S_hat = (S_hat+S_hat')/2; % adjust for rounding errors    

lden = -n*T/2*log(pi) -n/2*sum(h) -n/2*(sum(log(Hyper.VA)) +ldet(K_A))...
    +Hyper.nu0/2*ldet(Hyper.S0) -(Hyper.nu0+T)/2*ldet(S_hat) ...
    +mgammaln(n,(Hyper.nu0+T)/2) -mgammaln(n,Hyper.nu0/2);
end