tmpyhat0 = zeros(nsims,2*n+1);  %% [point forecasts, prelike, joint den]
tmpyhat1 = zeros(nsims,2*n+1); 
tmpyhat2 = zeros(nsims,2*n+1); 
tmpyhat4 = zeros(nsims,2*n+1); 

%% compute a few things
XX = X'*X;
KA = sparse(1:k,1:k,1./VA0) + XX;
Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + X'*shortYt);
Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + shortYt'*shortYt ...
    - Ahat'*KA*Ahat;
Shat = (Shat+Shat')/2; % adjust for rounding errors
nuhat = nu0+Tt;

for isim = 1:nsims + burnin
  
        %% sample Sig and A  
    Sig = iwishrnd(Shat,nuhat);
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig'; 
    
    if isim > burnin
        isave = isim - burnin;
        xtp1 = [1 reshape(shortYt(end:-1:end-p+1,:)',1,n*p)];
        if is_last_miss % if the lastest data are missing, do one more forecast horizon
            Ytp1 = xtp1*A + (CSig*randn(n,1))';
            xtp1 = [1 Ytp1 xtp1(2:end-n)];
        end            
        for tt=1:5
            EYtp1 = xtp1*A;
            dSig = diag(Sig)';
            if tt == 1
                tmpu = CSig\(data_tpk(1,:)-EYtp1)';
                lden_joint = -n/2*log(2*pi) -sum(diag(log(CSig))) -.5*(tmpu'*tmpu);
                lden = -.5*log(2*pi*dSig) - .5*(data_tpk(1,:)-EYtp1).^2./dSig;                
                tmpyhat0(isave,:) = [EYtp1 lden lden_joint]; 
            elseif tt == 2 && t<=T-tt
                tmpu = CSig\(data_tpk(2,:)-EYtp1)';
                lden_joint = -n/2*log(2*pi) -sum(diag(log(CSig))) -.5*(tmpu'*tmpu);
                lden = -.5*log(2*pi*dSig) - .5*(data_tpk(2,:)-EYtp1).^2./dSig;                
                tmpyhat1(isave,:) = [EYtp1 lden lden_joint]; 
            elseif tt == 3 && t<=T-tt
                tmpu = CSig\(data_tpk(3,:)-EYtp1)';
                lden_joint = -n/2*log(2*pi) -sum(diag(log(CSig))) -.5*(tmpu'*tmpu);
                lden = -.5*log(2*pi*dSig) - .5*(data_tpk(3,:)-EYtp1).^2./dSig;                
                tmpyhat2(isave,:) = [EYtp1 lden lden_joint]; 
            elseif tt == 5 && t<=T-tt
                tmpu = CSig\(data_tpk(5,:)-EYtp1)';
                lden_joint = -n/2*log(2*pi) -sum(diag(log(CSig))) -.5*(tmpu'*tmpu);
                lden = -.5*log(2*pi*dSig) - .5*(data_tpk(5,:)-EYtp1).^2./dSig;                
                tmpyhat4(isave,:) = [EYtp1 lden lden_joint];
            end
            Ytp1 = EYtp1 + (CSig*randn(n,1))';
            xtp1 = [1 Ytp1 xtp1(2:end-n)];
        end        
        
    end
end

