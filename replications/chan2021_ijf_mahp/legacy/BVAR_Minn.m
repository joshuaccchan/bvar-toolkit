% This function estimates a large Bayesian VAR with a data-based Minnesota 
% prior
%
% See:
% Chan, J.C.C. (2021). Minnesota-Type Adaptive Hierarchical Priors for 
% Large Bayesian VARs, International Journal of Forecasting, forthcoming

    % initialize the Markov chain
kappa = [.04,.04,1,100];
h0 = log(sig2);
h = repmat(h0',T,1);
Sigh = Sh0;
beta = zeros(k_beta,1);
alp = zeros(k_alp,1);

store_kappa = zeros(nsim,2);
store_beta = zeros(nsim,k_beta);
store_alp = zeros(nsim,k_alp);
store_Sigh = zeros(nsim,n);
store_h = zeros(nsim,T,n);

start_time = clock;

    % MCMC starts here
randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
disp('Starting MCMC for BVAR-Minn.... ');

for isim = 1:nsim + burnin
        % sample alp and beta    
    [Valp,Vbeta] = getVtheta(idx_kappa1,idx_kappa2,kappa,C,sig2);
    count_alp = 0;
    U = zeros(T,n);
    for ii = 1:n
        yi = Y(:,ii);
        ki = n*p+ii;
        Xi = [Z -Y(:,1:ii-1)];
        
        iVthetai = sparse(1:ki,1:ki,1./[Vbeta((ii-1)*(n*p+1)+1:ii*(n*p+1));...
            Valp(count_alp+1:count_alp+ii-1)]);
        XiiSighi = Xi'*sparse(1:T,1:T,exp(-h(:,ii)));
        Kthetai = iVthetai + XiiSighi*Xi;
        CKthetai = chol(Kthetai,'lower');
        thetai_hat = (CKthetai')\(CKthetai\(XiiSighi*yi));
        thetai = thetai_hat + CKthetai'\randn(ki,1);
        U(:,ii) = yi - Xi*thetai;
        
        beta((ii-1)*(n*p+1)+1:ii*(n*p+1)) =  thetai(1:n*p+1);
        alp(count_alp+1:count_alp+ii-1) =  thetai(n*p+2:end);
        count_alp = count_alp + ii -1;
    end
    
        % sample h
    for i=1:n
        Ystar = log(U(:,i).^2 + .0001);
        h(:,i) = SVRW(Ystar,h(:,i),Sigh(i),h0(i));
    end 
    
        % sample kappa1 and kappa2
    tmpc1 = sum(beta(idx_kappa1).^2./C(idx_kappa1));
    tmpc2 = sum(beta(idx_kappa2).^2./C(idx_kappa2));
    kappa(1) = gigrnd(c01(1)-n*p/2,2*c01(2),tmpc1,1);
    kappa(2) = gigrnd(c02(1)-(n-1)*n*p/2,2*c02(2),tmpc2,1);    
    
        % sample h0
    Kh0 = sparse(1:n,1:n,1./Sigh + 1./Vh);
    h0_hat = Kh0\(ah./Vh + h(1,:)'./Sigh);
    h0 = h0_hat + chol(Kh0,'lower')'\randn(n,1);
    
        % sample Sigh
    e = h - [h0';h(1:T-1,:)];
    Sigh = 1./gamrnd(nuh0+T/2,1./(Sh0 + sum(e.^2)'/2));    
    
    if isim > burnin        
        isave = isim-burnin;        
        store_kappa(isave,:) = kappa(1:2);
        store_Sigh(isave,:) = Sigh';
        store_beta(isave,:) =  beta';
        store_alp(isave,:) =  alp';
        store_h(isave,:,:) =  h;
    end
    
    if (mod(isim, 5000) == 0)
        disp([num2str(isim) ' loops... ' ])
    end     
end

disp( ['MCMC takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
disp(' ' );

h_hat = squeeze(mean(store_h));
alp_hat = mean(store_alp)';
beta_hat = mean(store_beta)';
kappa_hat = mean(store_kappa)';

figure;
subplot(1,2,1); 
histogram(store_kappa(:,1),50); 
xlabel('Histogram of posterior draws of \kappa_1'); box off;
subplot(1,2,2); 
histogram(store_kappa(:,2),50); 
xlabel('Histogram of posterior draws of \kappa_2'); box off;
set(gcf,'position',[300,300,800,400]);


