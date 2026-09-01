% This function computes the forecasts from a large Bayesian SVAR with the
% Minnesota-type Horseshoe prior, using CCCM algorithm

% construct a few useful things
tmpyhat = zeros(nsim,2*n+3,12);  %% [point forecasts, prelike, joint den]
k_alp = n*(n-1)/2;       % dimension of alp
k_beta = n^2*p + n;      % dimension of beta
k = k_beta/n;
sig2 = get_resid_var(Y0,Yt);
[C,idx_kappa1,idx_kappa2] = get_C(n,p,sig2);

% initialization for storeage
store_alp = zeros(nsim,k_alp);
store_beta = zeros(nsim,k_beta);
store_h_T = zeros(nsim,n);  % store only h_T
store_hpara = zeros(nsim,3*n);
store_kappa = zeros(nsim,2);
count_phi = zeros(n,1);

A_id = nonzeros(tril(reshape(1:n^2,n,n),-1)');
A = eye(n);

% initialize the Markov chain
sig2 = 1./gamrnd(Hyper.nuh,1./Hyper.Sh);
phi = Hyper.phi0; phi = min(Hyper.phi0 + sqrt(Hyper.Vphi).*randn(n,1),.99);
XX = Xt'*Xt; B = ((XX + .01*speye(k))\(Xt'*Yt))';
XB = Xt*B'; U = Yt-XB; Sig_hat = U'*U/Tt;
mu = zeros(n,1);
for ii=1:n
    s2i = U(:,ii).^2; mu(ii) = mean(log(s2i));
end
h = repmat(log(diag(Sig_hat))',Tt,1);
z_psi1 = 1./gamrnd(.5,1,n*p,1);
z_psi2 = 1./gamrnd(.5,1,(n-1)*n*p,1);
z_kappa = 1./gamrnd(.5,1,2,1);
psi_kappa1 = 1./gamrnd(.5,z_psi1);
psi_kappa2 = 1./gamrnd(.5,z_psi2);
Psi = ones(k*n,1);
Psi(idx_kappa1) = psi_kappa1; Psi(idx_kappa2) = psi_kappa2;

% MCMC starts here
for isim = 1:nsim + burnin
    % sample B
    tmpdV = getVbeta(idx_kappa1,idx_kappa2,kappa,C.*Psi,sig2);
    Ytilde = Yt*sparse(A');
    for ii = 1:n
        tmpXBA = XB(:,[1:ii-1 ii+1:n])*sparse(A(:,[1:ii-1 ii+1:n])');
        Zi = Ytilde(:,ii:n) - tmpXBA(:,ii:n);
        zi = reshape(Zi',Tt*(n-ii+1),1);
        Xi = repmat(Xt,n-ii+1,1);
        tmp1 = exp(-h(:,ii:n)).*repmat(A(ii:n,ii)',Tt,1);
        tmp2 = tmp1.*repmat(A(ii:n,ii)',Tt,1);
        iVbetai = sparse(1:k,1:k,1./tmpdV((ii-1)*k+1:ii*k));
        betai0 = Hyper.beta0((ii-1)*k+1:ii*k);
        Kbetai = iVbetai + Xi'*sparse(1:(n-ii+1)*Tt,1:(n-ii+1)*Tt,tmp2(:))*Xi;
        CKbetai = chol(Kbetai,'lower');
        betai_hat = (CKbetai')\(CKbetai\(iVbetai*betai0 ...
            + Xi'*sparse(1:(n-ii+1)*Tt,1:(n-ii+1)*Tt,tmp1(:))*Zi(:)));
        
        betai = betai_hat + CKbetai'\randn(k,1);
        B(ii,:) = betai;
        XB(:,ii) = Xt*betai;
    end
    beta = reshape(B',k_beta,1);
    
    % sample alp
    E = Yt - XB;
    count_alp = 0;
    for ii=2:n
        X_alpi = -E(:,1:ii-1);
        iD = sparse(1:Tt,1:Tt,exp(-h(:,ii)));
        iValpi = sparse(1:ii-1,1:ii-1,1./Hyper.Valp(count_alp+1:count_alp+ii-1));
        Kalpi = iValpi + X_alpi'*iD*X_alpi;
        alpi_hat = Kalpi\(X_alpi'*iD*E(:,ii));
        alpi = alpi_hat + chol(Kalpi,'lower')'\randn(ii-1,1);
        alp(count_alp+1:count_alp+ii-1) = alpi;
        count_alp = count_alp + ii-1;
    end
    A(A_id) = alp;
    
    % sample h
    AE = E*sparse(A'); 	
    for ii=1:n
        ystar = log(AE(:,ii).^2 + .0001);
        h(:,ii) = sample_SV(ystar,h(:,ii),0,phi(ii),sig2(ii));
    end
    
    % sample mu, phi and sig2
    [mu,phi,sig2,flag_phi] = sample_SVpara(h,mu,phi,Hyper);
    count_phi = count_phi + flag_phi;
    
    % sample psi
    tmpv1 = 1./z_psi1 + beta(idx_kappa1).^2./(2*C(idx_kappa1)*kappa(1));
    tmpv2 = 1./z_psi2 + beta(idx_kappa2).^2./(2*C(idx_kappa2)*kappa(2));
    psi_kappa1 = 1./gamrnd(1,1./tmpv1);
    psi_kappa2 = 1./gamrnd(1,1./tmpv2);
    Psi(idx_kappa1) = psi_kappa1;
    Psi(idx_kappa2) = psi_kappa2;
    
    % sample z_psi
    z_psi1 = 1./gamrnd(1,1./(1+1./psi_kappa1));
    z_psi2 = 1./gamrnd(1,1./(1+1./psi_kappa2));
    
    % sample kappa1 and kappa2
    tmpc1 = 1/z_kappa(1) + sum(beta(idx_kappa1).^2./(2*psi_kappa1.*C(idx_kappa1)));
    tmpc2 = 1/z_kappa(2) + sum(beta(idx_kappa2).^2./(2*psi_kappa2.*C(idx_kappa2)));
    kappa(1) = 1./gamrnd((n*p+1)/2,1./tmpc1);
    kappa(2) = 1./gamrnd(((n-1)*n*p+1)/2,1./tmpc2);
    
    % sample z_kappa
    z_kappa = 1./gamrnd(1,1./(1+1./kappa(1:2)));
    
    if isim > burnin
        isave = isim-burnin;
        store_beta(isave,:) = reshape(B',k_beta,1);
        store_alp(isave,:) = alp';
        store_h_T(isave,:) =  h(end,:)';
        store_hpara(isave,:) = [mu',phi',sig2'];
        store_kappa(isave,:) = kappa(1:2);
    end
end
kappa_hat = mean(store_kappa)';
kappaCI = quantile(store_kappa,[0.25 .975]);

insim=51;
% compute forecasts
for isim = 1:nsim
    tmpyhat_isim = zeros(insim,2*n+3,12);  %% [point forecasts, prelike, joint den]  
    for ijsim = 1:insim
        alp = store_alp(isim,:)';
        beta = store_beta(isim,:)';
        h_Tp1 = store_h_T(isim,:)';
        Sigh = store_hpara(isim,2*n+1:end)';
        phih = store_hpara(isim,n+1:2*n)';
        muh = store_hpara(isim,1:n)';
        
        % trasnform the parameters into reduced-form
        sqrtSigh = sqrt(Sigh);
        A(A_id) = alp;
        h_Tp1 = muh + phih.*(h_Tp1-muh) + sqrtSigh.*randn(n,1);
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
            
            h_Tp1 = muh + phih.*(h_Tp1-muh) + sqrtSigh.*randn(n,1);
            S = (A\sparse(1:n,1:n,exp(h_Tp1)))/A';
            CS = chol(S,'lower');
        end
    end
    tmpyhat(isim,:,:) = [squeeze(median(tmpyhat_isim(:,1:n,:))); squeeze(mean(tmpyhat_isim(:,n+1:end,:)))];   
end