% This script estimates the VAR-SV model in reduced-form with stationary
% autoregressive SV and computes the associated marginal likelihood
%
% See:
% Chan, J.C.C. (2022). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, forthcoming.

    % initalize storage
store_alp = zeros(nsim,k_alp);
store_beta = zeros(nsim,k_beta);
store_h = zeros(nsim,T,n);
store_hpara = zeros(nsim,3*n);
store_kappa = zeros(nsim,3);
count_phi = zeros(n,1);

    % initialize the Markov chain
mu = log(sig2_hat);
sig2 = 1./gamrnd(Hyper.nuh,1./Hyper.Sh);
phi = Hyper.phi0;
A = zeros(k,n);
h = zeros(T,n);
XX = X'*X;
for ii=1:n
    iValpi = sparse(1:k,1:k,1./Hyper.Valp((ii-1)*k+1:ii*k));
    A(:,ii) = (XX + iValpi)\(X'*Y(:,ii));
    s2i = (Y(:,ii) - X*A(:,ii)).^2;
    mu(ii) = mean(log(s2i));
    h(:,ii) = getARh_approx1N(s2i,mu(ii),phi(ii),sig2(ii));    
end
alp = reshape(A,k_alp,1);
beta = zeros(k_beta,1);
B0_id = nonzeros(tril(reshape(1:n^2,n,n),-1)');
B0 = eye(n);
[C_alp,idx_kappa1,idx_kappa2] = get_C(n,p,sig2_hat);
C_beta = Hyper.Vbeta/kappa4; % the constant part of Vbeta

    % MCMC starts here
randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
start_time = clock;
for isim = 1:nsim + burnin
        % sample alp
    [~,Hyper.Valp] = prior_Minn(p,kappa1,kappa2,kappa3,Y0,Y);    
    U = zeros(T,n);
    sqrt_exph = exp(h/2);
    for ii = 1:n        
        A(:,ii) = 0;
        Lambda = vec(sqrt_exph(:,ii:n));
        yi = vec((Y-X*A)*B0(ii:n,:)')./Lambda;    
        Wi = kron(B0(ii:n,ii),X)./Lambda;        
        iValpi = sparse(1:k,1:k,1./Hyper.Valp((ii-1)*k+1:ii*k));
        alpi0 = Hyper.alp0((ii-1)*k+1:ii*k);
        Kalpi = iValpi + Wi'*Wi;
        CKalpi = chol(Kalpi,'lower');        
        alpi_hat = (CKalpi')\(CKalpi\(iValpi*alpi0 + Wi'*yi));        
        alpi = alpi_hat + CKalpi'\randn(k,1);
        A(:,ii) = alpi;        
    end
    alp = reshape(A,k_alp,1);
    
        % sample beta
    [~,Hyper.Vbeta] = prior_B0(Y0,Y,kappa4);
    E = Y - X*A;
    count_beta = 0;
    for ii=2:n
        X_betai = -E(:,1:ii-1);   
        iD = sparse(1:T,1:T,exp(-h(:,ii)));
        iVbetai = sparse(1:ii-1,1:ii-1,1./Hyper.Vbeta(count_beta+1:count_beta+ii-1));         
        Kbetai = iVbetai + X_betai'*iD*X_betai;
        betai_hat = Kbetai\(X_betai'*iD*E(:,ii));
        betai = betai_hat + chol(Kbetai,'lower')'\randn(ii-1,1);
        beta(count_beta+1:count_beta+ii-1) = betai;        
        count_beta = count_beta + ii-1;
    end
    B0(B0_id) = beta;    
    
        % sample h
    B0E= E*B0';
    for ii=1:n
        ystar = log(B0E(:,ii).^2 + .0001);
        h(:,ii) = sample_SV(ystar,h(:,ii),mu(ii),phi(ii),sig2(ii)); 
    end    
    
        % sample mu, phi and sig2
    [mu,phi,sig2,flag_phi] = sample_SVpara(h,mu,phi,Hyper);
    count_phi = count_phi + flag_phi;   
    
        % sample kappa1, kappa2 and kappa4
    tmpc4 = sum(beta.^2./C_beta);
    kappa4 = gigrnd(Hyper.c0(3,1)-n*(n-1)/4,2*Hyper.c0(3,2),tmpc4,1);
    if is_kappafixed
        % nothing to do here
    elseif is_kappasym
        tmpc1 = sum(alp(idx_kappa1).^2./C_alp(idx_kappa1)) + ...
            sum(alp(idx_kappa2).^2./C_alp(idx_kappa2));
        kappa1 = gigrnd(Hyper.c0(1,1)-n^2*p/2,2*Hyper.c0(1,2),tmpc1,1);
        kappa2 = kappa1;
    else
        tmpc1 = sum(alp(idx_kappa1).^2./C_alp(idx_kappa1));
        tmpc2 = sum(alp(idx_kappa2).^2./C_alp(idx_kappa2));
        kappa1 = gigrnd(Hyper.c0(1,1)-n*p/2,2*Hyper.c0(1,2),tmpc1,1);
        kappa2 = gigrnd(Hyper.c0(2,1)-(n-1)*n*p/2,2*Hyper.c0(2,2),tmpc2,1);         
    end

    if ( mod(isim, 10000) ==0 )
        disp(  [ num2str(isim) ' loops... ' ] )
    end 
    
    if isim > burnin        
        isave = isim-burnin;         
        store_alp(isave,:) = reshape(A,1,k_alp);
        store_beta(isave,:) = beta';
        store_h(isave,:,:) = h;
        store_hpara(isave,:) = [mu',phi',sig2'];        
        store_kappa(isave,:) = [kappa1,kappa2,kappa4];
    end
end
disp( ['Estimating ' model_name ' takes '  num2str( etime( clock, start_time) ) ' seconds' ] );

alp_mean = mean(store_alp)';
A_mean = reshape(alp_mean,k,n)';
beta_mean = mean(store_beta)';
h_mean = squeeze(mean(store_h));
hpara_mean = mean(store_hpara)';
kappa_mean = mean(store_kappa)';
    
    % compute the marginal likelihood
if cp_ml
    disp(['Computing the marginal likelihood of ' model_name '... ']);
    start_time = clock;
    [lml,lmlstd] = ml_var_arsv_redu(X,Y,Y0,M,Hyper,flag_marg,store_h,...
        store_beta,store_hpara,store_kappa,is_kappafixed,is_kappasym);
    disp( ['ML computation takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
end
disp(' ' );
