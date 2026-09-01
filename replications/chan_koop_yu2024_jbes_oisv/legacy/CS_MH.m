% This is the main program for estimation of an CS-SV model
  
k_alp = n*(n-1)/2;       % dimension of alp
k_beta = n^2*p + n;      % dimension of beta
k = k_beta/n;

    % priors
sig2 = get_resid_var(Y0,Y);
kappa = [.1,.1,1,100]; % kappa(1): own lag; kappa(2): other lag; kappa(4): intercepts; kappa(3): alpha
[C,idx_kappa1,idx_kappa2] = get_C(n,p,sig2);
Hyper.beta0 = zeros(k_beta,1); 
%[Hyper.alp0,Hyper.Valp] = prior_alp(Y0,Y,kappa4);
Hyper.alp0 = zeros(n*(n-1)/2,1); Hyper.Valp = 1*ones(n*(n-1)/2,1); 
Hyper.mu0 = zeros(n,1); Hyper.Vmu = 100*ones(n,1);
Hyper.nuh = 3*ones(n,1); Hyper.Sh = .05*(Hyper.nuh-1);
Hyper.phi0 = .95*ones(n,1); Hyper.Vphi = .05^2*ones(n,1);

    % initalize storage
store_alp = zeros(nsim,k_alp);
store_beta = zeros(nsim,k_beta);
store_h = zeros(nsim,T,n);
store_hpara = zeros(nsim,3*n);
store_kappa = zeros(nsim,2);
count_phi = zeros(n,1);

A_id = nonzeros(tril(reshape(1:n^2,n,n),-1)');
A = eye(n);

    % initialize the Markov chain
sig2 = 1./gamrnd(Hyper.nuh,1./Hyper.Sh);
phi = Hyper.phi0; phi = min(Hyper.phi0 + sqrt(Hyper.Vphi).*randn(n,1),.99);
XX = X'*X; B = ((XX + .01*speye(k))\(X'*Y))'; 
XB = X*B'; U = Y-XB; Sig_hat = U'*U/T; 
mu = zeros(n,1);
for ii=1:n            
    s2i = U(:,ii).^2; mu(ii) = mean(log(s2i));   
end
h = repmat(log(diag(Sig_hat))',T,1); 
%beta = reshape(B',k_beta,1);
z_psi1 = 1./gamrnd(.5,1,n*p,1);
z_psi2 = 1./gamrnd(.5,1,(n-1)*n*p,1);
z_kappa = 1./gamrnd(.5,1,2,1);
psi_kappa1 = 1./gamrnd(.5,z_psi1);
psi_kappa2 = 1./gamrnd(.5,z_psi2);
Psi = ones(k*n,1);
Psi(idx_kappa1) = psi_kappa1; Psi(idx_kappa2) = psi_kappa2;

    % MCMC starts here
%randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
start_time = clock;
for isim = 1:nsim + burnin
        % sample B    
	tmpdV = getVbeta(idx_kappa1,idx_kappa2,kappa,C.*Psi,sig2);
    Ytilde = Y*sparse(A');   
    for ii = 1:n
        tmpXBA = XB(:,[1:ii-1 ii+1:n])*sparse(A(:,[1:ii-1 ii+1:n])');
        Zi = Ytilde(:,ii:n) - tmpXBA(:,ii:n);         
        zi = reshape(Zi',T*(n-ii+1),1);
        Xi = repmat(X,n-ii+1,1);
        tmp1 = exp(-h(:,ii:n)).*repmat(A(ii:n,ii)',T,1);
        tmp2 = tmp1.*repmat(A(ii:n,ii)',T,1);
        iVbetai = sparse(1:k,1:k,1./tmpdV((ii-1)*k+1:ii*k));
        betai0 = Hyper.beta0((ii-1)*k+1:ii*k);
        Kbetai = iVbetai + Xi'*sparse(1:(n-ii+1)*T,1:(n-ii+1)*T,tmp2(:))*Xi;            
        CKbetai = chol(Kbetai,'lower');
        betai_hat = (CKbetai')\(CKbetai\(iVbetai*betai0 ...
            + Xi'*sparse(1:(n-ii+1)*T,1:(n-ii+1)*T,tmp1(:))*Zi(:))); 

        betai = betai_hat + CKbetai'\randn(k,1);
        B(ii,:) = betai;
        XB(:,ii) = X*betai;
    end  
    beta = reshape(B',k_beta,1);
	
        % sample alp
    E = Y - XB;
    count_alp = 0;
    for ii=2:n
        X_alpi = -E(:,1:ii-1);   
        iD = sparse(1:T,1:T,exp(-h(:,ii)));
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
        h(:,ii) = sample_SV(ystar,h(:,ii),mu(ii),phi(ii),sig2(ii)); 
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
        store_h(isave,:,:) = h;
        store_hpara(isave,:) = [mu',phi',sig2']; 		
		store_kappa(isave,:) = kappa(1:2);
    end

    if (mod(isim, 2000) ==0)
        disp([num2str(isim) ' loops... ' ])
    end 	
end
disp( ['Estimation takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
disp(' ' );

beta_mean = mean(store_beta)'; Beta_mean = reshape(beta_mean,k,n)';
alp_mean = mean(store_alp)'; alp_median = median(store_alp)';
h_mean = squeeze(mean(store_h));
hpara_mean = mean(store_hpara)';
kappa_mean = mean(store_kappa)'; kappa_median = median(store_kappa)';


