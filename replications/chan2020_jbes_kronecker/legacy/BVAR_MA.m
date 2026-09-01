% This script estimates the BVAR-MA model
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

%% prior
psi0 = 0; Vpsi = 1;
lpri_psi = @(x) -.5*(x-psi0)^2/Vpsi -1e10*(x<-.99 || x>.99);
S0 = eye(n); nu0 = n+3;
construct_prior_A;

%% construct X
X = zeros(T,n*p); 
for i=1:p
    X(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end
X = [ones(T,1) X];
    
%% initialize for storage
store_Sig = zeros(n,n); 
store_A = zeros(k,n);
store_psi = zeros(nsims,1); 
ngrid = 300;
psigrid = linspace(-.1,.5,ngrid)';
store_ppsi = zeros(ngrid,1);
countpsi = 0;

%% initialize the chain
psi = .1;
Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
psihat = psi;

options = optimset('Display', 'off', 'LargeScale','off') ;

%% MCMC starts here
randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
disp('Starting MCMC for BVAR-MA.... ');
start_time = clock;

for isim = 1:nsims + burnin
  
        %% sample Sig and A    
    Xtld = Hpsi\X;    
    Ytld = Hpsi\shortY;
    iO = sparse(1:T,1:T,[1/(1+psi^2) ones(1,T-1)]);
    XtldiO = Xtld'*iO;
    KA = sparse(1:k,1:k,1./VA0) + XtldiO*Xtld;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XtldiO*Ytld);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + Ytld'*iO*Ytld ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    Sig = iwishrnd(Shat,nu0+T);
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig';   
    
    %% sample psi
    U = shortY - X*A;
    lp_psi = @(x) llike_MA(x,U,Sig) + lpri_psi(x);
    if (mod(isim,100)==0) || isim == 1; %% get the Hessian every 100 iterations
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
        alpMH =  lp_psi(psic) - lp_psi(psi) + ...
            -.5*(psi-psihat)^2*Kpsic + .5*(psic-psihat)^2*Kpsic;
    else
        alpMH = -inf;
    end
    if alpMH > log(rand)
        psi = psic;
        countpsi = countpsi + 1;
        Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
    end   
    
    if isim > burnin
        isave = isim - burnin; 
        store_A = store_A + A;
        store_Sig = store_Sig + Sig;        
        store_psi(isave,:) = psi;
        
%         % Uncomment this block to compute the posterior density of psi
%         % It would add a couple of minutes to the estimation
%         tmpden = zeros(ngrid,1);
%         for ii=1:ngrid
%             tmpden(ii) = lp_psi(psigrid(ii)); 
%         end
%         tmpden = exp(tmpden-max(tmpden));
%         tmpden = tmpden/(sum(tmpden)*(psigrid(2)-psigrid(1)));
%         store_ppsi = store_ppsi + tmpden;
    end
    
    if ( mod(isim, 5000) ==0 )
        disp(  [ num2str(isim) ' loops... ' ] )
    end 
    
end

disp( ['MCMC takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
disp(' ' );

A_mean = store_A/nsims;
Sig_mean = store_Sig/nsims;
psi_mean = mean(store_psi)';

if cp_ml
    ml_BVAR_MA;
end

figure;
colormap('hsv');
imagesc(A_mean);
colorbar;
box off;
title('Heat map of the VAR coefficients');    

figure; 
hist(store_psi,50);
title('Histogram of posterior draws of psi');
box off;

% figure; 
% plot(psigrid,ppsi_mean); box off;
% title('Posterior density of psi');
% box off;

