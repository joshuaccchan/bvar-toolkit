% This script estimates the marginal likelihood of the BVAR-CSV-MA model
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

disp('Computing the marginal likelihood of BVAR-CSV-MA... ');

llike = intlike_BVAR_CSV_MA(shortY,X,A_mean,Sig_mean,theta_mean(1),...
    theta_mean(2),theta_mean(3),5000);
ngrid = 300;

c_rho = 1/(normcdf(1,rho0,sqrt(Vrho))-normcdf(-1,rho0,sqrt(Vrho)));
c_psi = 1/(normcdf(1,psi0,sqrt(Vpsi))-normcdf(-1,psi0,sqrt(Vpsi)));

lpri = -.5*log(2*pi*Vpsi) + log(c_psi) -.5*(theta_mean(1)-psi0)^2/Vpsi ...
    -.5*log(2*pi*Vrho) + log(c_rho) -.5*(theta_mean(2)-rho0)^2/Vrho ...
    + nuh0*log(Sh0) - gammaln(nuh0) - (nuh0+1)*log(theta_mean(3)) - Sh0/theta_mean(3) ...
    + lniwpdf(A_mean,Sig_mean,A0,sparse(1:k,1:k,1./VA0),nu0,S0);

%% evaluate the posterior density
store_lpost = zeros(nsims,2); % [log density of A Sig, log density of sigh2]
for isim = 1:nsims
    h = store_h(isim,:)';
    psi = store_theta(isim,1);
    rho = store_theta(isim,2);
      
        %% compute the conditional density of Sig and A  
    Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
    Xtld = Hpsi\X;    
    Ytld = Hpsi\shortY;
    iO_hpsi = sparse(1:T,1:T,[1/(1+psi^2)*exp(-h(1)); exp(-h(2:end))]);
    XiO = Xtld'*iO_hpsi;
    KA = sparse(1:k,1:k,1./VA0) + XiO*Xtld;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiO*Ytld);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + Ytld'*iO_hpsi*Ytld ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    lden_ASig = lniwpdf(A_mean,Sig_mean,Ahat,KA,nu0+T,Shat);    
  
        %% compute the conditional density of sigh2
    eh = [h(1)*sqrt(1-rho^2);  h(2:end)-rho*h(1:end-1)];    
    lden_sigh2 = linvgammpdf(theta_mean(3),nuh0+T/2,Sh0 + sum(eh.^2)/2);
    
    store_lpost(isim,:) = [lden_ASig lden_sigh2];   
end
tmpmax = max(store_lpost);
lpost = log(mean(exp(store_lpost-repmat(tmpmax,nsims,1)))) + tmpmax;

nsims2 = 1000;
store_lpost = zeros(nsims2,2); % [log density of psi, log density of rho]
rhogrid = sort([theta_mean(2); linspace(-.999,.999,ngrid)']);
rhoidx = find(rhogrid==theta_mean(2));
psigrid = sort([theta_mean(1); linspace(-.99,.99,ngrid)']);
psiidx = find(psigrid==theta_mean(1));

U = shortY - X*A_mean;
CSig = chol(Sig_mean,'lower');
sigh2 = theta_mean(3);
ht = h_mean;
for isim = 1:nsims2    
    Utld = Hpsi\U;
    tmp = (Utld/CSig');
    s2 = sum(tmp.^2,2);
    s2(1) = s2(1)/(1+psi^2);
    h = sample_h(s2,rho,sigh2,h,n);
    
        %% sample rho
    Krho = 1/Vrho + sum(h(1:T-1).^2)/sigh2;
    rhohat = Krho\(rho0/Vrho + h(1:T-1)'*h(2:T)/sigh2);
    rhoc = rhohat + sqrt(Krho)'\randn;
    grho = @(x) -.5*log(sigh2./(1-x.^2))-.5*(1-x.^2)/sigh2*h(1)^2;
    if abs(rhoc)<.9999
        alpMH = exp(grho(rhoc)-grho(rho));
        if alpMH>rand
            rho = rhoc;            
            Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
        end
    end
    
    %% sample psi
    lp_psi = @(x) llike_CSV_MA(x,U,Sig_mean,h) + lpri_psi(x);    
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
        alpMH =  lp_psi(psic) - lp_psi(psi) + ...
            -.5*(psi-psihat)^2*Kpsic + .5*(psic-psihat)^2*Kpsic;
    else
        alpMH = -inf;
    end
    if alpMH > log(rand)
        psi = psic;        
        Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
    end  
  
        %% compute the conditional density of rho
    tmpden = grho(rhogrid) + -.5*Krho*(rhogrid-rhohat).^2;
    tmpden = exp(tmpden-max(tmpden));    
    tmpden = tmpden/(sum(tmpden)*(rhogrid(2)-rhogrid(1)));
    den_rho = tmpden(rhoidx);     
    
        %% compute the conditional density of psi
    tmpden = zeros(ngrid+1,1);
    for ii=1:ngrid+1
        tmpden(ii) = lp_psi(psigrid(ii)); 
    end
    tmpden = exp(tmpden-max(tmpden));
    tmpden = tmpden/(sum(tmpden)*(psigrid(2)-psigrid(1)));
    den_psi = tmpden(psiidx);

    store_lpost(isim,:) = [den_rho den_psi];   
end
lpost(3:4) = log(mean(store_lpost));

ML = llike + lpri - sum(lpost);
fprintf('log marginal likelihood of BVAR-CSV-MA: %.1f \n', ML); 
