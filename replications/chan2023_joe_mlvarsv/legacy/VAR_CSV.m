% This script estimates the VAR-CSV model and computes the associated 
% marginal likelihood
%
% See:
% Chan, J.C.C. (2022). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, forthcoming.

    % initialize for storage
store_Sig = zeros(nsim,n*(n+1)/2); 
store_a = zeros(nsim,k*n);
store_h = zeros(nsim,T);
store_hpara = zeros(nsim,2);
store_kappa = zeros(nsim,1);

    % initialize the chain
sig2 = 1./gamrnd(Hyper.nuh,1./Hyper.Sh);
phi = Hyper.phi0;
[Hyper.A0,Hyper.VA] = prior_NCP(p,kappa,kappa3,Y0,Y);
C = Hyper.VA/kappa; % the constant part of VA
K_A = sparse(1:k,1:k,1./Hyper.VA) + X'*X;
A = K_A\(sparse(1:k,1:k,Hyper.VA)\Hyper.A0 + X'*Y);
Sig = (Hyper.S0 + Hyper.A0'*sparse(1:k,1:k,1./Hyper.VA)*Hyper.A0 ...
        + Y'*Y - A'*K_A*A)/T;
U = Y - X*A;
tmp = (U/chol(Sig,'lower')');
s2 = sum(tmp.^2,2);
h = sample_CSV(s2,phi,sig2,zeros(T,1),n,true);
Hphi = speye(T) - phi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
count_h = 0; count_phi = 0;

    % MCMC starts here
randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
start_time = clock;
for isim = 1:nsim + burnin  
        % sample Sig and A  
    [~,Hyper.VA] = prior_NCP(p,kappa,kappa3,Y0,Y);
    iOh = sparse(1:T,1:T,exp(-h));
    XiOh = X'*iOh;
    K_A = sparse(1:k,1:k,1./Hyper.VA) + XiOh*X;
    A_hat = K_A\(sparse(1:k,1:k,Hyper.VA)\Hyper.A0 + XiOh*Y);
    S_hat = Hyper.S0 + Hyper.A0'*sparse(1:k,1:k,1./Hyper.VA)*Hyper.A0 ...
        + Y'*iOh*Y - A_hat'*K_A*A_hat;
    S_hat = (S_hat+S_hat')/2; % adjust for rounding errors
    Sig = iwishrnd(S_hat,Hyper.nu0+T);
    CSig = chol(Sig,'lower');
    A = A_hat + (chol(K_A,'lower')'\randn(k,n))*CSig'; 
    
        % sample h
    U = Y - X*A;
    tmp = (U/CSig');
    s2 = sum(tmp.^2,2);
    if isim < 20
        h = sample_CSV(s2,phi,sig2,h,n,true);
        flag = 0;
    else
        [h,flag] = sample_CSV(s2,phi,sig2,h,n);
    end
    count_h = count_h + flag;
    
        % sample phi and sig2
    [~,phi,sig2,flag_phi] = sample_SVpara(h,0,phi,Hyper);
    count_phi = count_phi + flag_phi;
    Hphi = speye(T) - phi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
    
        % sample kappa  
    if ~is_kappafixed
        Q = diag((A-Hyper.A0)*(Sig\(A-Hyper.A0)'));
        tmpc = sum(Q(2:end)./C(2:end));
        kappa = gigrnd(Hyper.c0(1)-n^2*p/2,2*Hyper.c0(2),tmpc,1); % k=np+1
    end
    
    if isim > burnin
        isave = isim - burnin; 
        store_a(isave,:,:) = A(:);
        store_Sig(isave,:) = vech(Sig);
        store_h(isave,:) = h';
        store_hpara(isave,:) = [phi sig2];
        store_kappa(isave,:) = kappa;
    end
    
    if ( mod(isim, 10000) ==0 )
        disp(  [ num2str(isim) ' loops... ' ] )
    end 
    
end
disp( ['Estimating ' model_name ' takes '  num2str( etime( clock, start_time) ) ' seconds' ] );

if cp_ml
    disp(['Computing the marginal likelihood of ' model_name '... ']);
    start_time = clock;
    [lml,lmlstd] = ml_var_csv(X,Y,Y0,M,Hyper,store_h,store_hpara,store_kappa,is_kappafixed);    
    disp( ['ML computation takes '  num2str(etime(clock, start_time)) ' seconds' ] );    
end
disp(' ' );

h_mean = mean(store_h)';
hpara_mean = mean(store_hpara)';
kappa_mean = mean(store_kappa);

CSV_std_mean = mean(exp(store_h/2))';

T_id = linspace(1961,2019.75,T)';
figure; 
plot(T_id, CSV_std_mean); box off;
title('Posterior mean of exp(h_t/2)');
box off; xlim([T_id(1)-1 T_id(end)+1]);
