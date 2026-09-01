% This script estimates the BVAR-t model
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

%% prior
S0 = eye(n); nu0 = n+3;
construct_prior_A;
nuub = 100; %% upperbound for nu

%% construct X
X = zeros(T,n*p); 
for i=1:p
    X(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end
X = [ones(T,1) X];
   
%% initialize for storage
store_Sig = zeros(n,n); 
store_A = zeros(k,n);
store_nu = zeros(nsims,1); 
store_lam = zeros(nsims,T);
nugrid = linspace(2,nuub,700)';
store_pnu = zeros(700,1);

%% initialize the chain
nu = 5;
lam = 1./gamrnd(nu/2,2/nu,T,1);
countnu = 0;

%% MCMC starts here
randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
disp('Starting MCMC for BVAR-t.... ');
start_time = clock;

for isim = 1:nsims + burnin
  
        %% sample Sig and A    
    iOm = sparse(1:T,1:T,1./lam);
    XiOm = X'*iOm;
    KA = sparse(1:k,1:k,1./VA0) + XiOm*X;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiOm*shortY);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + shortY'*iOm*shortY ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    Sig = iwishrnd(Shat,nu0+T);
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig'; 
    
        %% sample lam
    U = shortY - X*A;
    tmp = (U/CSig');
    s2 = sum(tmp.^2,2);
    lam = 1./gamrnd((n+nu)/2,2./(s2+nu));
    
        %% sample nu
    [nu,flag,fnu] = sample_nu(lam,nu,nuub);
    countnu = countnu + flag;

    if isim > burnin
        isave = isim - burnin; 
        store_A = store_A + A;
        store_Sig = store_Sig + Sig;        
        store_nu(isave,:) = nu;
        store_lam(isave,:) = lam';
        
        % compute the density of nu
        tmpden = fnu(nugrid);
        tmpden = exp(tmpden-max(tmpden));
        tmpden = tmpden/(sum(tmpden)*(nugrid(2)-nugrid(1)));
        store_pnu = store_pnu + tmpden; 
    end
    
    if ( mod(isim, 5000) ==0 )
        disp(  [ num2str(isim) ' loops... ' ] )
    end 
    
end

disp( ['MCMC takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
disp(' ' );

A_mean = store_A/nsims;
Sig_mean = store_Sig/nsims;
nu_mean = mean(store_nu)';
pnu_mean = store_pnu/nsims;

figure;
colormap('hsv');
imagesc(A_mean);
colorbar;
box off;
title('Heat map of the VAR coefficients');    

figure; 
plot(nugrid,pnu_mean); 
title('Posterior density of nu');
box off; xlim([0 30]);

if cp_ml
    ml_BVAR_t;
end