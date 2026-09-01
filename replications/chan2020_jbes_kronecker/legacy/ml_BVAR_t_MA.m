% This script estimates the marginal likelihood of the BVAR-t-MA model
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

disp('Computing the marginal likelihood of BVAR-t-MA... ');
psi_mean = theta_mean(1);
nu_mean = theta_mean(2);
ngrid = 300;

%% evaluate the log likelihood
Hpsi = speye(T) + psi_mean(1)*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
CSig = chol(Sig_mean,'lower');
Utld = Hpsi\(shortY - X*A_mean);
tmp = (Utld/CSig');
s2 = sum(tmp.^2,2);
s2(1) = s2(1)/(1+psi_mean^2);
llike = -T*n/2*log(nu_mean*pi) - n/2*log(1+psi_mean^2) ...
    + T*(gammaln((nu_mean+n)/2) - gammaln(nu_mean/2)) ...
    - T*sum(log(diag(CSig))) - (nu_mean+n)/2*sum(log(1+s2/nu_mean));

c_psi = 1/(normcdf(1,psi0,sqrt(Vpsi))-normcdf(-1,psi0,sqrt(Vpsi)));
lpri = lniwpdf(A_mean,Sig_mean,A0,sparse(1:k,1:k,1./VA0),nu0,S0) ...
    + log(1/(nuub-2)) ...
    -.5*log(2*pi*Vpsi) + log(c_psi) -.5*(psi_mean(1)-psi0)^2/Vpsi;

%% evaluate the posterior density
store_lpost = zeros(nsims,2); % [log density of A Sig, density of nu]
nugrid = sort([nu_mean; linspace(2,nuub,ngrid)']);
nuidx = find(nugrid==nu_mean);
for isim = 1:nsims
    lam = store_lam(isim,:)';
    psi1 = store_theta(isim,1);
    Hpsi = speye(T) + psi1*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
  
        %% compute the conditional density of Sig and A    
    Xtld = Hpsi\X;    
    Ytld = Hpsi\shortY;
    iO_lam = sparse(1:T,1:T,1./lam);
    XiO = Xtld'*iO_lam;
    KA = sparse(1:k,1:k,1./VA0) + XiO*Xtld;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiO*Ytld);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + Ytld'*iO_lam*Ytld ...
        - Ahat'*KA*Ahat;  
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    lden_ASig = lniwpdf(A_mean,Sig_mean,Ahat,KA,nu0+T,Shat);    
  
        %% sample nu
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
lpost = zeros(3,1);
lpost(1) = log(mean(exp(store_lpost(:,1)-tmpmax))) + tmpmax;
lpost(2) = log(mean(store_lpost(:,2)));

store_lpost = zeros(nsims,1); % density of psi
psigrid = sort([theta_mean(1); linspace(-.99,.99,ngrid)']);
psiidx = find(psigrid==theta_mean(1));

U = shortY - X*A_mean;
CSig = chol(Sig_mean,'lower');
for isim = 1:nsims 
    
    %% sample lam
    Utld = Hpsi\U;
    tmp = Utld/CSig';
    s2 = sum(tmp.^2,2);
    s2(1) = s2(1)/(1+psi1^2);
    lam = 1./gamrnd((n+nu_mean)/2,2./(s2+nu_mean));
    
    %% sample psi
    lp_psi = @(x) llike_CSV_MA(x,U,Sig_mean,log(lam)) + lpri_psi(x); 
    if (mod(isim,100)==0) || isim == 1; %% get the Hessian every 100 iterations
        [psihat,fval,exitflag,output,grad,hess] ...
            = fminunc(@(x)-lp_psi(x),psihat,options); 
        [tmpCpsi flag] = chol(hess,'lower');
        if flag == 0
            Kpsic = hess;
        else 
            Kpsic = 1/.05^2;
        end
    else
        psihat = fminbnd(@(x)-lp_psi(x),-.99,.99);
    end    
    psic = psihat + 1/sqrt(Kpsic)*randn;
    if abs(psic)<.99
        alpMH =  lp_psi(psic) - lp_psi(psi1) + ...
            -.5*(psi1-psihat)^2*Kpsic + .5*(psic-psihat)^2*Kpsic;
    else
        alpMH = -inf;
    end
    if alpMH > log(rand)
        psi1 = psic;        
        Hpsi = speye(T) + psi1*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
    end    
    
        %% compute the conditional density of psi
    tmpden = zeros(ngrid+1,1);
    for ii=1:ngrid+1
        tmpden(ii) = lp_psi(psigrid(ii)); 
    end
    tmpden = exp(tmpden-tmpden(psiidx));
    tmpden = tmpden/(sum(tmpden)*(psigrid(2)-psigrid(1)));
    den_psi = tmpden(psiidx);

    store_lpost(isim,:) = den_psi; 
end
lpost(3) = log(mean(store_lpost));
ML = llike + lpri - sum(lpost);
fprintf('log marginal likelihood of BVAR-t-MA: %.1f \n', ML); 
