tmpyhat0 = zeros(nsims,2*n+1);  %% [point forecasts, prelike, joint den]
tmpyhat1 = zeros(nsims,2*n+1); 
tmpyhat2 = zeros(nsims,2*n+1); 
tmpyhat4 = zeros(nsims,2*n+1); 

    % initialize
nu = 5;
lam = 1./gamrnd(nu/2,2/nu,Tt,1);

for isim = 1:nsims + burnin
  
        % sample Sig and A
    iOm = sparse(1:Tt,1:Tt,1./lam);
    XiOm = X'*iOm;
    KA = sparse(1:k,1:k,1./VA0) + XiOm*X;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiOm*shortYt);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + shortYt'*iOm*shortYt ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    Sig = iwishrnd(Shat,nu0+Tt);
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig'; 
    
        % sample lam
    U = shortYt - X*A;
    tmp = (U/CSig');
    s2 = sum(tmp.^2,2);
    lam = 1./gamrnd((n+nu)/2,2./(s2+nu));    
   
       % sample nu
    nu = sample_nu(lam,nu,nuub);   
       
    if isim > burnin
        isave = isim - burnin;
        xtp1 = [1 reshape(shortYt(end:-1:end-p+1,:)',1,n*p)];        
        if is_last_miss % if the lastest data are missing, do one more forecast horizon
            EYtp1 = xtp1*A;
            Ytp1 = EYtp1 + (CSig*randn(n,1))'/sqrt(gamrnd(nu/2,2/nu));
            xtp1 = [1 Ytp1 xtp1(2:end-n)];            
        end        
        for tt = 1:5
            EYtp1 = xtp1*A;
            dSig = diag(Sig)';
            ct = gammaln((nu+1)/2) - gammaln(nu/2) - .5*log(nu*pi*dSig);
            ct_joint = gammaln((nu+n)/2) - gammaln(nu/2) - n/2*log(nu*pi); 
            if tt == 1
                tmpu = CSig\(data_tpk(1,:)-EYtp1)';
                lden_joint = ct_joint - sum(diag(log(CSig)))...
                    -(nu+n)/2*log(1 + (tmpu'*tmpu)/nu);               
                lden = ct - (nu+1)/2*log(1 + (data_tpk(1,:)-EYtp1).^2./dSig/nu);   
                tmpyhat0(isave,:) = [EYtp1 lden lden_joint]; 
            elseif tt == 2 && t<=T-tt
                tmpu = CSig\(data_tpk(2,:)-EYtp1)';
                lden_joint = ct_joint - sum(diag(log(CSig))) ...
                    -(nu+n)/2*log(1 + (tmpu'*tmpu)/nu);               
                lden = ct - (nu+1)/2*log(1 + (data_tpk(2,:)-EYtp1).^2./dSig/nu);   
                tmpyhat1(isave,:) = [EYtp1 lden lden_joint]; 
            elseif tt == 3 && t<=T-tt
                tmpu = CSig\(data_tpk(3,:)-EYtp1)';
                lden_joint = ct_joint - sum(diag(log(CSig))) ...
                    -(nu+n)/2*log(1 + (tmpu'*tmpu)/nu);               
                lden = ct - (nu+1)/2*log(1 + (data_tpk(3,:)-EYtp1).^2./dSig/nu);   
                tmpyhat2(isave,:) = [EYtp1 lden lden_joint]; 
            elseif tt == 5 && t<=T-tt
                tmpu = CSig\(data_tpk(5,:)-EYtp1)';
                lden_joint = ct_joint - sum(diag(log(CSig))) ...
                    -(nu+n)/2*log(1 + (tmpu'*tmpu)/nu);               
                lden = ct - (nu+1)/2*log(1 + (data_tpk(5,:)-EYtp1).^2./dSig/nu);   
                tmpyhat4(isave,:) = [EYtp1 lden lden_joint];
            end              
            Ytp1 = EYtp1 + (CSig*randn(n,1))'/sqrt(gamrnd(nu/2,2/nu));
            xtp1 = [1 Ytp1 xtp1(2:end-n)];        
        end        
        
    end
end