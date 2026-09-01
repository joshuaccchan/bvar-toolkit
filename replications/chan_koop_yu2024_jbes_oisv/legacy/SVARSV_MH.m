% This is the main program for estimation of an OI-SV model

k = 1+n*p;
k_alpha = n*k;
sig2 = get_resid_var(Y0,Y);
[C,idx_kappa1,idx_kappa2] = get_C(n,p,sig2);

    % priors
Hyper.nuh = 3*ones(n,1); Hyper.Sh = .05*ones(n,1).*(Hyper.nuh-1);
Hyper.phi0 = .95*ones(n,1); Hyper.Vphi = .05^2*ones(n,1);
Hyper.B0 = eye(n); Hyper.VB0 = 1*ones(n); 

    % initialization for storeage
store_kappa = zeros(nsim,2);
store_alpha = zeros(nsim,k_alpha);
store_B0 = zeros(nsim,n^2);
store_h = zeros(nsim,T,n);
store_hpara = zeros(nsim,2*n);
count_phi = zeros(n,1);

start_time = clock;

%% MCMC starts here
randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
disp('Starting MCMC for BSVAR.... ');

    % initialize the Markov chain
kappa = [.1,.1,NaN,100]; % kappa(1): own lag; kappa(2): other lag; kappa(4): intercepts
A = (X'*X + .01*speye(k))\(X'*Y);
U = Y-X*A;
alpha = A(:);
Sig_hat = U'*U/T;
B0 = diag(1./sqrt(diag(Sig_hat)));
h = repmat(log(diag(Sig_hat))',T,1);
sig2 = 1./gamrnd(Hyper.nuh,1./Hyper.Sh);
phi = min(Hyper.phi0 + sqrt(Hyper.Vphi).*randn(n,1),.99);
E = zeros(T,n);
z_psi1 = 1./gamrnd(.5,1,n*p,1);
z_psi2 = 1./gamrnd(.5,1,(n-1)*n*p,1);
z_kappa = 1./gamrnd(.5,1,2,1);
psi_kappa1 = 1./gamrnd(.5,z_psi1);
psi_kappa2 = 1./gamrnd(.5,z_psi2);
Psi = ones(k*n,1);
Psi(idx_kappa1) = psi_kappa1; Psi(idx_kappa2) = psi_kappa2;

for isim = 1:nsim + burnin    
        % sammple B0
    U = Y-X*A;
    for ii=1:n
        EiOhi = U'*sparse(1:T,1:T,exp(-h(:,ii)));       
        Kbi = sparse(1:n,1:n,1./Hyper.VB0(ii,:)) + EiOhi*U; 
        mui = Kbi\(Hyper.B0(ii,:)./Hyper.VB0(ii,:))'; 
        Ci = chol(Kbi,'lower')/sqrt(T);
        Gam_mi = B0([1:ii-1 ii+1:end],:)';
        Gam_miperp = null(Gam_mi');
        
        V = zeros(n,n); zeta = zeros(n,1);        
        for jj=1:n
            if jj==1                
                v1 = Ci\Gam_miperp; v1 = v1/norm(v1);
                V = [v1 null(v1')];
                zetaj_hat = mui'*(Ci*v1);
                zeta(1) = anormrnd(zetaj_hat,1/T);
            else
                zetaj_hat = mui'*(Ci*V(:,jj));
                zeta(jj) = zetaj_hat + 1/sqrt(T)*randn;
            end
        end
        phii = (Ci')\sum(V.*repmat(zeta',n,1),2);
        % B0(ii,:) = phii;
        B0(ii,:) = phii*sign(phii(ii)); % fix the sign of the i-th element to be positive
    end   
    
        % sample alpha    
    tmpdV = getVbeta(idx_kappa1,idx_kappa2,kappa,C.*Psi,sig2);
    Lambda = vec(exp(h/2));
    for ii=1:n
        A(:,ii) = 0;
        yi = vec((Y-X*A)*B0')./Lambda;
        Wi = kron(B0(:,ii),X)./Lambda;       
        iValphai = sparse(1:k,1:k,1./tmpdV((ii-1)*k+1:ii*k));
        Kalphai = iValphai + Wi'*Wi; 
        CKalphai = chol(Kalphai,'lower');
        alphai_hat = (CKalphai')\(CKalphai\(Wi'*yi));        
        alphai = alphai_hat + CKalphai'\randn(k,1);
        A(:,ii) = alphai;         
    end  
    alpha = A(:);

        % sample h
    E = (Y - X*A)*B0';
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
        store_h(isave,:,:) =  h;
    end
    
    if (mod(isim, 2000) ==0)
        disp([num2str(isim) ' loops... ' ])
    end     
end

disp( ['MCMC takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
disp(' ' );

h_mean = squeeze(mean(store_h));
B0_mean = mean(store_B0)'; B0_median = median(store_B0)';
alpha_mean = mean(store_alpha)';
hpara_mean = mean(store_hpara)';
kappa_mean = mean(store_kappa)'; kappa_median = median(store_kappa)';
A_mean = reshape(alpha_mean,k,n)';

