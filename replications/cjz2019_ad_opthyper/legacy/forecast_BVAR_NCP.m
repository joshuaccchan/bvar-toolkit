tmpyhat1 = zeros(nsims,2*n+1);  %% [point forecasts, prelike, joint den]
tmpyhat4 = zeros(nsims,2*n+1); 

%% construct the natural conjugate prior
if model == 1
    [A0,VA,nu0,S0] = prior_NCP(p,kappa_0,Y0,Yt);
elseif model == 2
    kappa_con_opt = fmincon(@(kappa)ThreeKappa(kappa,kappa_0(1,4:5),Yt,Y0,Tt,n,k,p),...
    kappa_con_opt,[],[],[],[],zeros(3,1),[Inf;10;Inf],[],options);
    [A0,VA,nu0,S0] = prior_NCP(p,[kappa_con_opt kappa_0(1,4:5)],Y0,Yt);
elseif model == 3
    kappa_opt = fmincon(@(kappa)FiveKappa(kappa,Yt,Y0,Tt,n,k,p),...
    kappa_opt,[],[],[],[],zeros(5,1),[Inf;10;Inf;Inf;Inf],[],options);
    [A0,VA,nu0,S0] = prior_NCP(p,kappa_opt,Y0,Yt);
end
iVA = VA\speye(k);

%% compute a few things
tmpY = [Y0(end-p+1:end,:); Yt];
Z = zeros(Tt,n*p);
for i=1:p
    Z(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end
Z = [ones(Tt,1) Z];

ZZ = Z'*Z;
KA = iVA + ZZ;
Ahat = KA\(iVA*A0 + Z'*Yt);
Shat = S0 + A0'*iVA*A0 + Yt'*Yt - Ahat'*KA*Ahat;
Shat = (Shat+Shat')/2; % adjust for rounding errors
nuhat = nu0+Tt;

for isim = 1:nsims + burnin
  
        % sample Sig and A  
    Sig = iwishrnd(Shat,nuhat);
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig'; 
    
    if isim > burnin
        isave = isim - burnin;
        xtp1 = [1 reshape(Yt(end:-1:end-p+1,:)',1,n*p)];
        for tt=1:4
            EYtp1 = xtp1*A;
            dSig = diag(Sig)';
            if tt == 1
                tmpu = CSig\(Y(t+1,:)-EYtp1)';
                lden_joint = -n/2*log(2*pi) -sum(diag(log(CSig))) -.5*(tmpu'*tmpu);
                lden = -.5*log(2*pi*dSig) - .5*(Y(t+1,:)-EYtp1).^2./dSig;                
                tmpyhat1(isave,:) = [EYtp1 lden lden_joint]; 
            elseif tt == 4 && t<=T-tt
                tmpu = CSig\(Y(t+4,:)-EYtp1)';
                lden_joint = -n/2*log(2*pi) -sum(diag(log(CSig))) -.5*(tmpu'*tmpu);
                lden = -.5*log(2*pi*dSig) - .5*(Y(t+4,:)-EYtp1).^2./dSig;                
                tmpyhat4(isave,:) = [EYtp1 lden lden_joint];          
            end
            Ytp1 = EYtp1 + (CSig*randn(n,1))';
            xtp1 = [1 Ytp1 xtp1(2:end-n)];
        end        
        
    end
end

