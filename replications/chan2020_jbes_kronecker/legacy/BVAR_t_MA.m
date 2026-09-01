% This script estimates the BVAR-t-MA model
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

%% prior
psi0 = 0; Vpsi = 1;
lpri_psi = @(x) -.5*(x-psi0)^2/Vpsi -1e10*(x<-.99 || x>.99);
S0 = eye(n); nu0 = n+3;
nuub = 100; %% upperbound for nu
construct_prior_A;

%% construct X
X = zeros(T,n*p); 
for i=1:p
    X(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end
X = [ones(T,1) X];

%% initialize for storeage
store_Sig = zeros(n,n); 
store_A = zeros(k,n);
store_lam = zeros(nsims,T);
store_theta = zeros(nsims,2);

options = optimset('Display', 'off', 'LargeScale','off') ;
countnu = 0; countpsi = 0;

ngrid = 300;
psigrid = linspace(-.1,.5,ngrid)';
store_ppsi = zeros(ngrid,1);
nugrid = linspace(2,50,ngrid)';
store_pnu = zeros(ngrid,1);

%% initialize the chain
nu = 5;
lam = 1./gamrnd(nu/2,2/nu,T,1);
psi1 = -.1;
Hpsi = speye(T) + psi1*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
psihat = psi1;

%% MCMC starts here
randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
disp('Starting MCMC for BVAR-t-MA.... ');
start_time = clock;

for isim = 1:nsims + burnin
  
        %% sample Sig and A   
    Xtld = Hpsi\X;    
    Ytld = Hpsi\shortY;
    iO_lam = sparse(1:T,1:T,1./lam);
    XiO = Xtld'*iO_lam;
    KA = sparse(1:k,1:k,1./VA0) + XiO*Xtld;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiO*Ytld);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + Ytld'*iO_lam*Ytld ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    Sig = iwishrnd(Shat,nu0+T);    
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig';     

      %% sample psi1
    U = shortY - X*A;
    lp_psi = @(x) llike_CSV_MA(x,U,Sig,log(lam)) + lpri_psi(x);    
    if (mod(isim,100)==0) || isim == 1 %% get the Hessian every 100 iterations
        [psihat,fval,exitflag,output,grad,hess] ...
            = fminunc(@(x)-lp_psi(x),psihat,options); 
        [tmpCpsi flag] = chol(hess,'lower');
        if flag == 0
            Kpsic = hess;
        else 
            Kpsic = 1/.05^2;
        end
    else
        psihat = fminbnd(@(x)-lp_psi(x),-.99,.99);
    end    
    psic = psihat + 1/sqrt(Kpsic)*randn;
    if abs(psic)<.99
        alpMH =  lp_psi(psic) - lp_psi(psi1) + ...
            -.5*(psi1-psihat)^2*Kpsic + .5*(psic-psihat)^2*Kpsic;
    else
        alpMH = -inf;
    end
    if alpMH > log(rand)
        psi1 = psic;
        countpsi = countpsi + 1;
        Hpsi = speye(T) + psi1*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
    end     
    
        %% sample lam
    Utld = Hpsi\(shortY - X*A);
    tmp = (Utld/CSig');
    s2 = sum(tmp.^2,2);
    s2(1) = s2(1)/(1+psi1^2);
    lam = 1./gamrnd((n+nu)/2,2./(s2+nu));
    
        %% sample nu           
    [nu flag fnu] = sample_nu(lam,nu,nuub);
    countnu = countnu + flag;    
  
    if isim > burnin
        isave = isim - burnin; 
        store_A = store_A + A;
        store_Sig = store_Sig + Sig;
        store_lam(isave,:) = lam';
        store_theta(isave,:) = [psi1 nu];
        
%         % Uncomment this block to compute the posterior density of psi
%         % It would add a couple of minutes to the estimation
%         tmpden = zeros(ngrid,1);
%         for ii=1:ngrid
%             tmpden(ii) = lp_psi(psigrid(ii)); 
%         end
%         tmpden = exp(tmpden-max(tmpden));
%         tmpden = tmpden/(sum(tmpden)*(psigrid(2)-psigrid(1)));
%         store_ppsi = store_ppsi + tmpden;
%         
        % compute the posterior density of nu
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
theta_mean = mean(store_theta)';
ppsi_mean = store_ppsi/nsims;
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

figure; 
hist(store_theta(:,1),50);
title('Histogram of posterior draws of psi');
box off;

% figure; 
% plot(psigrid,ppsi_mean); box off;
% title('Posterior density of psi');
% box off;

if cp_ml
    ml_BVAR_t_MA;
end






