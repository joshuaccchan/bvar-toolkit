% This script produces forecasts from a 20-variable BVAR with a CSV
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
Sig = S0;

%% initialize
h = zeros(Tt,1);
rho = .8;
sigh2 = .1;
Hrho = speye(Tt) - rho*sparse(2:Tt,1:(Tt-1),ones(1,Tt-1),Tt,Tt);

for isim = 1:nsims + burnin
  
        % sample Sig and A    
    iOh = sparse(1:Tt,1:Tt,exp(-h));
    ZiOh = Z'*iOh;
    KA = sparse(1:k,1:k,1./VA0) + ZiOh*Z;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + ZiOh*shortYt);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + shortYt'*iOh*shortYt ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    Sig = iwishrnd(Shat,nu0+Tt);    
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig'; 
    
        % sample h
    U = shortYt - Z*A;
    tmp = (U/CSig');
    s2_h = sum(tmp.^2,2);    
    h = sample_h(s2_h,rho,sigh2,h,n);
    
        % sample sigh2
    eh = [h(1)*sqrt(1-rho^2);  h(2:end)-rho*h(1:end-1)];    
    sigh2 = 1/gamrnd(nuh0+Tt/2,1/(Sh0 + sum(eh.^2)/2));

        % sample rho
    Krho = 1/Vrho + sum(h(1:Tt-1).^2)/sigh2;
    rhohat = Krho\(rho0/Vrho + h(1:Tt-1)'*h(2:Tt)/sigh2);
    rhoc = rhohat + sqrt(Krho)'\randn;
    grho = @(x) -.5*log(sigh2./(1-x.^2))-.5*(1-x.^2)/sigh2*h(1)^2;
    if abs(rhoc)<.999
        alpMH = exp(grho(rhoc)-grho(rho));
        if alpMH>rand
            rho = rhoc;            
            Hrho = speye(Tt) - rho*sparse(2:Tt,1:(Tt-1),ones(1,Tt-1),Tt,Tt);
        end
    end    
       
    if isim > burnin
        isave = isim - burnin;
        xtp1 = [1 reshape(shortYt(end:-1:end-p+1,:)',1,n*p)];
        htp1 = h(end); 
        if is_last_miss % if the lastest data are missing, do one more forecast horizon
            htp1 = rho*htp1 + sqrt(sigh2)*randn;
            EYtp1 = xtp1*A;
            Ytp1 = EYtp1 + (exp(htp1/2)*CSig*randn(n,1))';            
            xtp1 = [1 Ytp1 xtp1(2:end-n)];
        end        
        for tt = 1:2
            htp1 = rho*htp1 + sqrt(sigh2)*randn;
            EYtp1 = xtp1*A;
            dSig = exp(htp1)*diag(Sig)';            
            if tt == 1
                tmpu = CSig\(data_tpk(1,:)-EYtp1)';
                lden_joint = -n/2*log(2*pi) -n/2*htp1 -sum(diag(log(CSig)))...
                    -.5*(tmpu'*tmpu)/exp(htp1);
                lden = -.5*log(2*pi*dSig) - .5*(data_tpk(1,:)-EYtp1).^2./dSig;                
                tmpyhat0(isave,:) = [EYtp1 lden lden_joint]; 
            elseif tt == 2 && t<=T-tt
                tmpu = CSig\(data_tpk(2,:)-EYtp1)';
                lden_joint = -n/2*log(2*pi) -n/2*htp1 -sum(diag(log(CSig)))...
                    -.5*(tmpu'*tmpu)/exp(htp1);
                lden = -.5*log(2*pi*dSig) - .5*(data_tpk(2,:)-EYtp1).^2./dSig;                
                tmpyhat1(isave,:) = [EYtp1 lden lden_joint];            
            end         
            Ytp1 = EYtp1 + (exp(htp1/2)*CSig*randn(n,1))';
            xtp1 = [1 Ytp1 xtp1(2:end-n)];
        end        
        
    end
end