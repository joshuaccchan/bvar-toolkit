% This script estimates the marginal likelihood of the BVAR-t model
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

disp('Computing the marginal likelihood of BVAR-t... ');

A_mean = store_A/nsims;
Sig_mean = store_Sig/nsims;
nu_mean = mean(store_nu)';

%% evaluate the log likelihood
CSig = chol(Sig_mean,'lower');
tmp = (shortY-X*A_mean)/CSig';
s2 = sum(tmp.^2,2);
llike = T*(gammaln((nu_mean+n)/2) - gammaln(nu_mean/2) - n/2*log(nu_mean*pi))...
    - T*sum(log(diag(CSig))) -(nu_mean+n)/2*sum(log(1+s2/nu_mean));
lpri = lniwpdf(A_mean,Sig_mean,A0,sparse(1:k,1:k,1./VA0),nu0,S0) + log(1/(nuub-2));

%% evaluate the posterior density
store_lpost = zeros(nsims,2); % [log density of A Sig, density of nu]
nugrid = sort([nu_mean; linspace(2,nuub,700)']);
nuidx = find(nugrid==nu_mean);

for isim = 1:nsims
    lam = store_lam(isim,:)';
  
        %% compute the conditional density of Sig and A    
    iOm = sparse(1:T,1:T,1./lam);
    XiOm = X'*iOm;
    KA = sparse(1:k,1:k,1./VA0) + XiOm*X;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiOm*shortY);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + shortY'*iOm*shortY ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    lden_ASig = lniwpdf(A_mean,Sig_mean,Ahat,KA,nu0+T,Shat);    
  
        %% compute the conditional density of nu
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

fprintf('log marginal likelihood of BVAR-t: %.1f \n', ML); 


