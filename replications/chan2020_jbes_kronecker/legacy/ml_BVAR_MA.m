% This script estimates the marginal likelihood of the BVAR-MA model
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

disp('Computing the marginal likelihood of BVAR-MA... ');

%% evaluate the log likelihood
Hpsi = speye(T) + psi_mean*sparse(2:T,1:(T-1),ones(1,T-1),T,T); 
CSig = chol(Sig_mean,'lower');
Utld = Hpsi\(shortY-X*A_mean);
tmp = (Utld/CSig');
s2 = sum(tmp.^2,2);
llike =  -T*n/2*log(2*pi) - T*sum(log(diag(CSig))) -n/2*log(1+psi_mean^2) ...
    -.5*(s2(1)/(1+psi^2) + sum(s2(2:end)));
c_psi = 1/(normcdf(1,psi0,sqrt(Vpsi))-normcdf(-1,psi0,sqrt(Vpsi)));
lpri = lniwpdf(A_mean,Sig_mean,A0,sparse(1:k,1:k,1./VA0),nu0,S0) ...
    -.5*log(2*pi*Vpsi) + log(c_psi) -.5*(psi_mean(1)-psi0)^2/Vpsi;

%% evaluate the posterior density
store_lpost = zeros(nsims,1); % [log density of A Sig]

for isim = 1:nsims
    psi = store_psi(isim,:)';
    Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T); 
  
        %% compute the conditional density of Sig and A    
    Xtld = Hpsi\X;    
    Ytld = Hpsi\shortY;
    iO = sparse(1:T,1:T,[1/(1+psi^2) ones(1,T-1)]);
    XtldiO = Xtld'*iO;
    KA = sparse(1:k,1:k,1./VA0) + XtldiO*Xtld;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XtldiO*Ytld);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + Ytld'*iO*Ytld ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    lden_ASig = lniwpdf(A_mean,Sig_mean,Ahat,KA,nu0+T,Shat);
     
    store_lpost(isim,:) = lden_ASig; 
end
tmpmax = max(store_lpost(:,1));
lpost = log(mean(exp(store_lpost(:,1)-tmpmax))) + tmpmax;

U = shortY - X*A_mean;
lp_psi = @(x) llike_MA(x,U,Sig_mean) + lpri_psi(x);
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
fprintf('log marginal likelihood of BVAR-MA: %.1f \n', ML);
