% This script estimates the marginal likelihood of the BVAR-CSV model
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

disp('Computing the marginal likelihood of BVAR-CSV... ');

llike = intlike_BVAR_CSV(shortY,X,A_mean,Sig_mean,theta_mean(1),theta_mean(2),1000);
c_rho = 1/(normcdf(1,rho0,sqrt(Vrho))-normcdf(-1,rho0,sqrt(Vrho)));
lpri = -.5*log(2*pi*Vrho) + log(c_rho) -.5*(theta_mean(1)-rho0)^2/Vrho ...
    + nuh0*log(Sh0) - gammaln(nuh0) - (nuh0+1)*log(theta_mean(2)) - Sh0/theta_mean(2) ...
    + lniwpdf(A_mean,Sig_mean,A0,sparse(1:k,1:k,1./VA0),nu0,S0);

%% evaluate the posterior density
store_lpost = zeros(nsims,2); % [log density of A Sig, log density of sigh2]

for isim = 1:nsims
    h = store_h(isim,:)';
    rho = store_theta(isim,1);
  
        %% compute the conditional density of Sig and A         
    iOh = sparse(1:T,1:T,exp(-h));
    XiOh = X'*iOh;
    KA = sparse(1:k,1:k,1./VA0) + XiOh*X;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiOh*shortY);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + shortY'*iOh*shortY ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors    
    lden_ASig = lniwpdf(A_mean,Sig_mean,Ahat,KA,nu0+T,Shat);    
  
        %% compute the conditional density of sigh2
    eh = [h(1)*sqrt(1-rho^2);  h(2:end)-rho*h(1:end-1)];    
    lden_sigh2 = linvgammpdf(theta_mean(2),nuh0+T/2,Sh0 + sum(eh.^2)/2);
    
    store_lpost(isim,:) = [lden_ASig lden_sigh2];   
end
tmpmax = max(store_lpost);
lpost = log(mean(exp(store_lpost-repmat(tmpmax,nsims,1)))) + tmpmax;

store_lpost = zeros(nsims,1); % [log density of rho]
rhogrid = sort([theta_mean(1); linspace(-.999,.999,700)']);
rhoidx = find(rhogrid==theta_mean(1));
U = shortY - X*A_mean;
CSig = chol(Sig_mean,'lower');
tmp = (U/CSig');
s2 = sum(tmp.^2,2);
sigh2 = theta_mean(2);
for isim = 1:nsims    
    HiSH = Hrho'*sparse(1:T,1:T,[(1-rho^2)/sigh2; 1/sigh2*ones(T-1,1)])*Hrho;
    errh = 1; ht = h_mean;
    while errh> 10^(-3);
        eht = exp(ht);
        sieht = s2./eht; 
        fh = -n/2 + .5*sieht;
        Gh = .5*sieht;
        Kh = HiSH + sparse(1:T,1:T,Gh);
        newht = Kh\(fh+Gh.*ht);
        errh = max(abs(newht-ht));
        ht = newht;          
    end
    CKh = chol(Kh,'lower');     
    % AR-step
    hstar = ht;    
    logc = -.5*hstar'*HiSH*hstar - n/2*sum(hstar) - .5*exp(-hstar)'*s2 + log(3);
    flag = 0;
    while flag == 0
        hc = ht + CKh'\randn(T,1);        
        alpARc = -.5*hc'*HiSH*hc - n/2*sum(hc) - .5*exp(-hc)'*s2 ...
            + .5*(hc-ht)'*Kh*(hc-ht) - logc;
        if alpARc > log(rand)
            flag = 1;
        end
    end        
    % MH-step    
    alpAR = -.5*h'*HiSH*h - n/2*sum(h) -.5*exp(-h)'*s2 ...
        + .5*(h-ht)'*Kh*(h-ht) - logc;
    if alpAR < 0
        alpMH = 1;
    elseif alpARc < 0
        alpMH = - alpAR;
    else
        alpMH = alpARc - alpAR;
    end    
    if alpMH > log(rand) || isim == 1
        h = hc;        
    end
    
        %% sample rho
    Krho = 1/Vrho + sum(h(1:T-1).^2)/sigh2;
    rhohat = Krho\(rho0/Vrho + h(1:T-1)'*h(2:T)/sigh2);
    rhoc = rhohat + sqrt(Krho)'\randn;
    grho = @(x) -.5*log(sigh2./(1-x.^2))-.5*(1-x.^2)/sigh2*h(1)^2;
    if abs(rhoc)<.9999
        alpMH = exp(grho(rhoc)-grho(rho));
        if alpMH>rand
            rho = rhoc;
            countrho = countrho+1;
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
lpost(3) = log(mean(store_lpost));
ML = llike + lpri - sum(lpost);
fprintf('log marginal likelihood of BVAR-CSV: %.1f \n', ML);
