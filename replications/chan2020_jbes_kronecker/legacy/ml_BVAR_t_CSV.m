% This script estimates the marginal likelihood of the BVAR-t-CSV model
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

disp('Computing the marginal likelihood of BVAR-t-CSV... ');

nu_mean = theta_mean(1);
rho_mean = theta_mean(2);
sigh2_mean = theta_mean(3);
llike = intlike_BVAR_t_CSV(shortY,X,A_mean,Sig_mean,rho_mean,sigh2_mean,...
    nu_mean,1000);
c_rho = 1/(normcdf(1,rho0,sqrt(Vrho))-normcdf(-1,rho0,sqrt(Vrho)));
lpri = log(1/(nuub-2)) ...
    -.5*log(2*pi*Vrho) + log(c_rho) -.5*(rho_mean-rho0)^2/Vrho ...
    + nuh0*log(Sh0) - gammaln(nuh0) - (nuh0+1)*log(sigh2_mean) - Sh0/sigh2_mean ...
    + lniwpdf(A_mean,Sig_mean,A0,sparse(1:k,1:k,1./VA0),nu0,S0);

%% evaluate the posterior density
store_lpost = zeros(nsims,3); % [log density of A Sig, log density of sigh2, log density of nu]
nugrid = sort([nu_mean; linspace(2,nuub,700)']);
nuidx = find(nugrid==nu_mean);

for isim = 1:nsims
    h = store_h(isim,:)';
    lam = store_lam(isim,:)';
    rho = store_theta(isim,2);
  
        %% compute the conditional density of Sig and A         
    iOm = sparse(1:T,1:T,exp(-h)./lam);
    XiOm = X'*iOm;
    KA = sparse(1:k,1:k,1./VA0) + XiOm*X;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiOm*shortY);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + shortY'*iOm*shortY ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    lden_ASig = lniwpdf(A_mean,Sig_mean,Ahat,KA,nu0+T,Shat);    
  
        %% compute the conditional density of sigh2
    eh = [h(1)*sqrt(1-rho^2);  h(2:end)-rho*h(1:end-1)];    
    lden_sigh2 = linvgammpdf(sigh2_mean,nuh0+T/2,Sh0 + sum(eh.^2)/2);
    
        %% compute the conditional density of nu
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

store_lpost = zeros(nsims,1); % [log density of rho]
rhogrid = sort([rho_mean; linspace(-.999,.999,700)']);
rhoidx = find(rhogrid==rho_mean);
U = shortY - X*A_mean;
CSig = chol(Sig_mean,'lower');
tmp = U/CSig';
s2 = sum(tmp.^2,2);
nu = nu_mean;
sigh2 = sigh2_mean;
for isim = 1:nsims    
    s2 = sum(tmp.^2,2)./lam;
    h = sample_h(s2,rho,sigh2,h,n);
    
        %% sample lam    
    s2 = sum(tmp.^2,2)./exp(h);
    lam = 1./gamrnd((n+nu)/2,2./(s2+nu));    
    
        %% sample rho
    Krho = 1/Vrho + sum(h(1:T-1).^2)/sigh2;
    rhohat = Krho\(rho0/Vrho + h(1:T-1)'*h(2:T)/sigh2);
    rhoc = rhohat + sqrt(Krho)'\randn;
    grho = @(x) -.5*log(sigh2./(1-x.^2))-.5*(1-x.^2)/sigh2*h(1)^2;
    if abs(rhoc)<.999
        alpMH = exp(grho(rhoc)-grho(rho));
        if alpMH>rand
            rho = rhoc;            
            Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
        end
    end 
  
        %% compute the conditional density of rho
    tmpden = grho(rhogrid) + -.5*Krho*(rhogrid-rhohat).^2;
    tmpden = exp(tmpden-max(tmpden));    
    tmpden = tmpden/(sum(tmpden)*(rhogrid(2)-rhogrid(1)));
    den_rho = tmpden(rhoidx);     
    store_lpost(isim,:) = den_rho;   
end
lpost(4) = log(mean(store_lpost));
ML = llike + lpri - sum(lpost);
fprintf('log marginal likelihood of BVAR-t-CSV: %.1f \n', ML); 