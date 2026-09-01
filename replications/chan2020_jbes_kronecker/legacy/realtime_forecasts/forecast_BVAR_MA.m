tmpyhat0 = zeros(nsims,2*n+1);  %% [point forecasts, prelike, joint den]
tmpyhat1 = zeros(nsims,2*n+1); 
tmpyhat2 = zeros(nsims,2*n+1); 
tmpyhat4 = zeros(nsims,2*n+1); 

%% initialize
psi = .1;
Hpsi = speye(Tt) + psi*sparse(2:Tt,1:(Tt-1),ones(1,Tt-1),Tt,Tt);
options = optimset('Display', 'off', 'LargeScale','off') ;
psihat = psi;
for isim = 1:nsims + burnin
  
        % sample Sig and A    
    Xtld = Hpsi\X;    
    Ytld = Hpsi\shortYt;
    iO = sparse(1:Tt,1:Tt,[1/(1+psi^2) ones(1,Tt-1)]);
    XtldiO = Xtld'*iO;
    KA = sparse(1:k,1:k,1./VA0) + XtldiO*Xtld;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XtldiO*Ytld);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + Ytld'*iO*Ytld ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    Sig = iwishrnd(Shat,nu0+Tt);
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig';   
    
        % sample psi
    U = shortYt - X*A;
    lp_psi = @(x) llike_MA(x,U,Sig) + lpri_psi(x);
    if (mod(isim,100)==0) || isim == 1 %% get the Hessian every 100 iterations
        [psihat,fval,exitflag,output,grad,hess] ...
            = fminunc(@(x)-lp_psi(x),psihat,options); 
        [tmpCpsi,flag] = chol(hess,'lower');
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
        Hpsi = speye(Tt) + psi*sparse(2:Tt,1:(Tt-1),ones(1,Tt-1),Tt,Tt);
    end   
       
    if isim > burnin
        isave = isim - burnin;
        xtp1 = [1 reshape(shortYt(end:-1:end-p+1,:)',1,n*p)];
        E = Hpsi\(shortYt - X*A);
        etp1 = E(end,:)';
        if is_last_miss % if the lastest data are missing, do one more forecast horizon            
            EYtp1 = xtp1*A + psi*etp1';
            etp1 = CSig*randn(n,1);
            Ytp1 = EYtp1 + etp1';            
            xtp1 = [1 Ytp1 xtp1(2:end-n)];
        end        
        for tt = 1:5
            EYtp1 = xtp1*A + psi*etp1';
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
            etp1 = CSig*randn(n,1);
            Ytp1 = EYtp1 + etp1';
            xtp1 = [1 Ytp1 xtp1(2:end-n)];
        end        
        
    end
end