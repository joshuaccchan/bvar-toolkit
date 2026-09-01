% This function computes the forecasts from a large Bayesian SVAR with the
% Minnesota-type Horseshoe prior, SVARSV

% construct a few useful things
tmpyhat = zeros(nsim,2*n+3,12);  %% [point forecasts, prelike, joint den]
k = 1+n*p;
k_alpha = n*k;
sig2 = get_resid_var(Y0,Yt);
[C,idx_kappa1,idx_kappa2] = get_C(n,p,sig2);

% initialization for storeage
store_kappa = zeros(nsim,2);
store_alpha = zeros(nsim,k_alpha);
store_B0 = zeros(nsim,n^2);
store_h_T = zeros(nsim,n);  % store only h_T
store_hpara = zeros(nsim,2*n);
count_phi = zeros(n,1);

% initialize the Markov chain
kappa = [.1,.1,NaN,100]; % kappa(1): own lag; kappa(2): other lag; kappa(4): intercepts
A = (Xt'*Xt + .01*speye(k))\(Xt'*Yt);
XA = Xt*A;
U = Yt-XA;
alpha = A(:);
Sig_hat = U'*U/Tt;
B0 = diag(1./sqrt(diag(Sig_hat)));
h = repmat(log(diag(Sig_hat))',Tt,1);
sig2 = 1./gamrnd(Hyper.nuh,1./Hyper.Sh);
phi = min(Hyper.phi0 + sqrt(Hyper.Vphi).*randn(n,1),.99);
E = zeros(Tt,n);
z_psi1 = 1./gamrnd(.5,1,n*p,1);
z_psi2 = 1./gamrnd(.5,1,(n-1)*n*p,1);
z_kappa = 1./gamrnd(.5,1,2,1);
psi_kappa1 = 1./gamrnd(.5,z_psi1);
psi_kappa2 = 1./gamrnd(.5,z_psi2);
Psi = ones(k*n,1);
Psi(idx_kappa1) = psi_kappa1; Psi(idx_kappa2) = psi_kappa2;

% MCMC starts here
for isim = 1:nsim + burnin
    % sammple B0
    U = Yt-Xt*A;
    for ii=1:n
        EiOhi = U'*sparse(1:Tt,1:Tt,exp(-h(:,ii)));
        Kbi = sparse(1:n,1:n,1./Hyper.VB0(ii,:)) + EiOhi*U; 
        mui = Kbi\(Hyper.B0(ii,:)./Hyper.VB0(ii,:))';
        Ci = chol(Kbi,'lower')/sqrt(Tt);
        Gam_mi = B0([1:ii-1 ii+1:end],:)';
        Gam_miperp = null(Gam_mi');
        
        V = zeros(n,n); zeta = zeros(n,1);
        for jj=1:n
            if jj==1
                v1 = Ci\Gam_miperp; v1 = v1/norm(v1);
                V = [v1 null(v1')];
                zetaj_hat = mui'*(Ci*v1);
                zeta(1) = anormrnd(zetaj_hat,1/Tt);
            else
                zetaj_hat = mui'*(Ci*V(:,jj));
                zeta(jj) = zetaj_hat + 1/sqrt(Tt)*randn;
            end
        end
        phii = (Ci')\sum(V.*repmat(zeta',n,1),2);
        % B0(ii,:) = phii;
        B0(ii,:) = phii*sign(phii(ii)); % fix the sign of the i-th element to be positive
    end
    
    % sample alpha
    tmpdV = getVbeta(idx_kappa1,idx_kappa2,kappa,C.*Psi,sig2);
    for ii=1:n
        zi = reshape(B0*(Yt-[XA(:,1:ii-1),sparse(Tt,1),XA(:,ii+1:end)])',Tt*n,1);
        Wi = kron(Xt,B0(:,ii));
        iValphai = sparse(1:k,1:k,1./tmpdV((ii-1)*k+1:ii*k));
        WiiD = Wi'*sparse(1:n*Tt,1:n*Tt,exp(-reshape(h',Tt*n,1)));
        Kalphai = iValphai + WiiD*Wi;
        CKalphai = chol(Kalphai,'lower');
        alphai_hat = (CKalphai')\(CKalphai\(WiiD*zi));
        alphai = alphai_hat + CKalphai'\randn(k,1);
        A(:,ii) = alphai;
        XA(:,ii) = Xt*alphai;
    end
    alpha = A(:);
    
    % sample h
    E = (Yt - XA)*B0';
    for ii=1:n
        ystar = log(E(:,ii).^2 + .0001);
        h(:,ii) = sample_SV(ystar,h(:,ii),0,phi(ii),sig2(ii));
    end
    
    % sample phi and sig2
    [phi,sig2,flag_phi] = sample_SV0para(h,phi,Hyper);
    count_phi = count_phi + flag_phi;
    
    % sample psi
    tmpv1 = 1./z_psi1 + alpha(idx_kappa1).^2./(2*C(idx_kappa1)*kappa(1));
    tmpv2 = 1./z_psi2 + alpha(idx_kappa2).^2./(2*C(idx_kappa2)*kappa(2));
    psi_kappa1 = 1./gamrnd(1,1./tmpv1);
    psi_kappa2 = 1./gamrnd(1,1./tmpv2);
    Psi(idx_kappa1) = psi_kappa1;
    Psi(idx_kappa2) = psi_kappa2;
    
    % sample z_psi
    z_psi1 = 1./gamrnd(1,1./(1+1./psi_kappa1));
    z_psi2 = 1./gamrnd(1,1./(1+1./psi_kappa2));
    
    % sample kappa1 and kappa2
    tmpc1 = 1/z_kappa(1) + sum(alpha(idx_kappa1).^2./(2*psi_kappa1.*C(idx_kappa1)));
    tmpc2 = 1/z_kappa(2) + sum(alpha(idx_kappa2).^2./(2*psi_kappa2.*C(idx_kappa2)));
    kappa(1) = 1./gamrnd((n*p+1)/2,1./tmpc1);
    kappa(2) = 1./gamrnd(((n-1)*n*p+1)/2,1./tmpc2);
    
    % sample z_kappa
    z_kappa = 1./gamrnd(1,1./(1+1./kappa(1:2)));
    
    if isim > burnin
        isave = isim-burnin;
        store_kappa(isave,:) = kappa(1:2);
        store_hpara(isave,:) = [phi',sig2'];
        store_alpha(isave,:) =  alpha';
        store_B0(isave,:) =  reshape(B0',n^2,1); % stacked by rows
        store_h_T(isave,:) =  h(end,:)';
    end
end
kappa_hat = mean(store_kappa)';
kappaCI = quantile(store_kappa,[0.25 .975]);


insim=51;
% compute forecasts
for isim = 1:nsim
    tmpyhat_isim = zeros(insim,2*n+3,12);  %% [point forecasts, prelike, joint den]     
    for ijsim = 1:insim
        alp = store_B0(isim,:)';
        beta = store_alpha(isim,:)';
        h_Tp1 = store_h_T(isim,:)';
        Sigh = store_hpara(isim,n+1:end)';
        phih = store_hpara(isim,1:n)';
        
        % trasnform the parameters into reduced-form
        sqrtSigh = sqrt(Sigh);
        A = reshape(alp,n,n)'; 
        h_Tp1 = phih.*h_Tp1 + sqrtSigh.*randn(n,1);
        S = (A\sparse(1:n,1:n,exp(h_Tp1)))/A';
        B = reshape(beta,n*p+1,n); 
        xtp1 = [1 reshape(Yt(end:-1:end-p+1,:)',1,n*p)];
        CS = chol(S,'lower');
        for tt=1:12
            EYtp1 = xtp1*B;
            dS = diag(S)';
			if t<=T-tt
			    tmpu = CS\(Y(t+tt,:)-EYtp1)';
                lden_joint = -n/2*log(2*pi) -sum(diag(log(CS))) -.5*(tmpu'*tmpu);
				lden_joint_first4 = -4/2*log(2*pi) -sum(diag(log(CS(1:4,1:4)))) -.5*(tmpu(1:4)'*tmpu(1:4));
				lden_joint_last4 = -4/2*log(2*pi) -sum(diag(log(CS(end-3:end,end-3:end)))) -.5*(tmpu(end-3:end)'*tmpu(end-3:end));
                lden = -.5*log(2*pi*dS) - .5*(Y(t+tt,:)-EYtp1).^2./dS;
                tmpyhat_isim(ijsim,:,tt) = [EYtp1 lden lden_joint lden_joint_first4 lden_joint_last4];
			end
            
            Ytp1 = EYtp1 + (CS*randn(n,1))';
            xtp1 = [1 Ytp1 xtp1(2:end-n)];
            
            h_Tp1 = phih.*h_Tp1 + sqrtSigh.*randn(n,1);
            S = (A\sparse(1:n,1:n,exp(h_Tp1)))/A';
            CS = chol(S,'lower');
        end
    end
    tmpyhat(isim,:,:) = [squeeze(median(tmpyhat_isim(:,1:n,:))); squeeze(mean(tmpyhat_isim(:,n+1:end,:)))];  
end