% This script produces forecasts from a 20-variable BVAR with a CSV + t +
% MA errors
%
% See:
% Chan, J.C.C. (2020). Large Bayesian Vector Autoregressions. In: P. Fuleky (Eds),
% Macroeconomic Forecasting in the Era of Big Data, 95-125, Springer, Cham

tmpyhat0 = zeros(nsims,2*n+1);  %% [point forecasts, prelike, joint den]
tmpyhat1 = zeros(nsims,2*n+1); 

%% construct the natural conjugate prior
[A0,VA0,nu0,S0] = prior_NC(p,c1,c2,Y0,shortYt);

tmpY = [Y0(end-p+1:end,:); shortYt];
Z = zeros(Tt,n*p);
for i=1:p
    Z(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end
Z = [ones(Tt,1) Z];

%% initialize
h = zeros(Tt,1);
nu = 5;
rho = .8;
sigh2 = .1;
Hrho = speye(Tt) - rho*sparse(2:Tt,1:(Tt-1),ones(1,Tt-1),Tt,Tt);
lam = 1./gamrnd(nu/2,2/nu,Tt,1);
psi = .1;
Hpsi = speye(Tt) + psi*sparse(2:Tt,1:(Tt-1),ones(1,Tt-1),Tt,Tt);
options = optimset('Display', 'off', 'LargeScale','off') ;
psihat = psi;

for isim = 1:nsims + burnin  
        % sample Sig and A   
    Ztld = Hpsi\Z;    
    Ytld = Hpsi\shortYt;
    iO_h_lam_psi = sparse(1:Tt,1:Tt,[1/(1+psi^2)*exp(-h(1));exp(-h(2:end))]./lam);
    ZiO = Ztld'*iO_h_lam_psi;
    KA = sparse(1:k,1:k,1./VA0) + ZiO*Ztld;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + ZiO*Ytld);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + Ytld'*iO_h_lam_psi*Ytld ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    Sig = iwishrnd(Shat,nu0+Tt);    
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig'; 
    
        % sample lam
    U = shortYt - Z*A;
    Utld = Hpsi\U;
    tmp = (Utld/CSig');
    s2_lam = sum(tmp.^2,2)./exp(h);
    lam = 1./gamrnd((n+nu)/2,2./(s2_lam+nu));   
    
        % sample h    
    s2_h = sum(tmp.^2,2)./lam;
    s2_h(1) = s2_h(1)/(1+psi^2);    
    h = sample_h(s2_h,rho,sigh2,h,n);    
    
        % sample sigh2
    eh = [h(1)*sqrt(1-rho^2);  h(2:end)-rho*h(1:end-1)];    
    sigh2 = 1/gamrnd(nuh0+Tt/2,1/(Sh0 + sum(eh.^2)/2));

        % sample rho
    Krho = 1/Vrho + sum(h(1:Tt-1).^2)/sigh2;
    rhohat = Krho\(rho0/Vrho + h(1:Tt-1)'*h(2:Tt)/sigh2);
    rhoc = rhohat + sqrt(Krho)'\randn;
    grho = @(x) -.5*log(sigh2./(1-x.^2))-.5*(1-x.^2)/sigh2*h(1)^2;
    if abs(rhoc)<.99
        alpMH = exp(grho(rhoc)-grho(rho));
        if alpMH>rand
            rho = rhoc;            
            Hrho = speye(Tt) - rho*sparse(2:Tt,1:(Tt-1),ones(1,Tt-1),Tt,Tt);
        end
    end 
    
        % sample psi
    U_psi = U./repmat(sqrt(lam),1,n);
    lp_psi = @(x) llike_CSV_MA(x,U_psi,Sig,h) + lpri_psi(x);    
    if (mod(isim,100)==0) || isim == 1 %% get the Hessian every 100 iterations
        [psihat,fval,exitflag,output,grad,hess] ...
            = fminunc(@(x)-lp_psi(x),psihat,options); 
        [tmpCpsi, flag] = chol(hess,'lower');
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
    
        % sample nu
    nu = sample_nu(lam,nu,nuub);    
       
    if isim > burnin
        isave = isim - burnin;
        xtp1 = [1 reshape(shortYt(end:-1:end-p+1,:)',1,n*p)];        
        E = Hpsi\(shortYt - Z*A);
        etp1 = E(end,:)';
        htp1 = h(end);        
        if is_last_miss % if the lastest data are missing, do one more forecast horizon
            htp1 = rho*htp1 + sqrt(sigh2)*randn;
            EYtp1 = xtp1*A + psi*etp1';
            etp1 = (exp(htp1/2)*CSig*randn(n,1))/sqrt(gamrnd(nu/2,2/nu));
            Ytp1 = EYtp1 + etp1';
            xtp1 = [1 Ytp1 xtp1(2:end-n)];
        end                
        for tt = 1:2
            htp1 = rho*htp1 + sqrt(sigh2)*randn;
            EYtp1 = xtp1*A + psi*etp1';
            dSig = exp(htp1)*diag(Sig)';
            ct = gammaln((nu+1)/2) - gammaln(nu/2) - .5*log(nu*pi*dSig);
            ct_joint = gammaln((nu+n)/2) - gammaln(nu/2) - n/2*log(nu*pi);
            if tt == 1
                tmpu = CSig\(data_tpk(1,:)-EYtp1)';
                lden_joint = ct_joint - sum(diag(log(CSig))) -n/2*htp1...
                    -(nu+n)/2*log(1 + (tmpu'*tmpu)/nu/exp(htp1));               
                lden = ct - (nu+1)/2*log(1 + (data_tpk(1,:)-EYtp1).^2./dSig/nu);   
                tmpyhat0(isave,:) = [EYtp1 lden lden_joint]; 
            elseif tt == 2 && t<=T-tt
                tmpu = CSig\(data_tpk(2,:)-EYtp1)';
                lden_joint = ct_joint - sum(diag(log(CSig))) -n/2*htp1...
                    -(nu+n)/2*log(1 + (tmpu'*tmpu)/nu/exp(htp1));
                lden = ct - (nu+1)/2*log(1 + (data_tpk(2,:)-EYtp1).^2./dSig/nu);   
                tmpyhat1(isave,:) = [EYtp1 lden lden_joint]; 
            end
            etp1 = (exp(htp1/2)*CSig*randn(n,1))/sqrt(gamrnd(nu/2,2/nu)); 
            Ytp1 = EYtp1 + etp1';
            xtp1 = [1 Ytp1 xtp1(2:end-n)];
        end        
    end
end