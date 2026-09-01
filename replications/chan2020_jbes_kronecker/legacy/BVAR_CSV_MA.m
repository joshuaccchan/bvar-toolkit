% This script estimates the BVAR-CSV-MA model
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

%% prior
psi0 = 0; Vpsi = 1;
lpri_psi = @(x) -.5*(x-psi0)^2/Vpsi -1e10*(x<-.99 || x>.99);
S0 = eye(n); nu0 = n+3;
nuh0 = 5; Sh0 = .01*(nuh0-1); 
rho0 = .9; Vrho = .2^2;
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
store_h = zeros(nsims,T);
store_theta = zeros(nsims,3);
options = optimset('Display', 'off', 'LargeScale','off') ;
counth = 0; countrho = 0; countpsi = 0;
ngrid = 300;
psigrid = linspace(-.1,.4,ngrid)';
store_ppsi = zeros(ngrid,1);

%% initialize the chain
h = zeros(T,1);
psi = -.1;
rho = .8;
sigh2 = .1;
Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
Hpsi = speye(T) + psi*sparse(2:T,1:(T-1),ones(1,T-1),T,T);

psihat = psi;

%% MCMC starts here
randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
disp('Starting MCMC for BVAR-CSV-MA.... ');
start_time = clock;

for isim = 1:nsims + burnin
  
        %% sample Sig and A   
    Xtld = Hpsi\X;    
    Ytld = Hpsi\shortY;
    iO_hpsi = sparse(1:T,1:T,[1/(1+psi^2)*exp(-h(1)); exp(-h(2:end))]);
    XiO = Xtld'*iO_hpsi;
    KA = sparse(1:k,1:k,1./VA0) + XiO*Xtld;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiO*Ytld);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + Ytld'*iO_hpsi*Ytld ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    Sig = iwishrnd(Shat,nu0+T);    
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig'; 
    
        %% sample h
    U = shortY - X*A;
    Utld = Hpsi\U;
    tmp = (Utld/CSig');
    s2 = sum(tmp.^2,2);
    s2(1) = s2(1)/(1+psi^2);
    [h flag] = sample_h(s2,rho,sigh2,h,n);
    counth = counth + flag;   
    
        %% sample sigh2
    eh = [h(1)*sqrt(1-rho^2);  h(2:end)-rho*h(1:end-1)];    
    sigh2 = 1/gamrnd(nuh0+T/2,1/(Sh0 + sum(eh.^2)/2));

        %% sample rho
    Krho = 1/Vrho + sum(h(1:T-1).^2)/sigh2;
    rhohat = Krho\(rho0/Vrho + h(1:T-1)'*h(2:T)/sigh2);
    rhoc = rhohat + sqrt(Krho)'\randn;
    grho = @(x) -.5*log(sigh2./(1-x.^2))-.5*(1-x.^2)/sigh2*h(1)^2;
    if abs(rhoc)<.9999
        alpMH = exp(grho(rhoc)-grho(rho));
        if alpMH>rand
            rho = rhoc;
            countrho = countrho+1;
            Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
        end
    end 
    
      %% sample psi
    U = shortY - X*A;
    lp_psi = @(x) llike_CSV_MA(x,U,Sig,h) + lpri_psi(x);    
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
        store_h(isave,:) = h';
        store_theta(isave,:) = [psi rho sigh2];
        
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
h_mean = mean(exp(store_h))';
theta_mean = mean(store_theta)';
ppsi_mean = store_ppsi/nsims;
CSV_std_mean = mean(exp(store_h/2))';

figure;
colormap('hsv');
imagesc(A_mean);
colorbar;
box off;
title('Heat map of the VAR coefficients');    

T_id = linspace(1959,2013.75,T)';
figure; 
plot(T_id, CSV_std_mean); box off;
title('Posterior mean of exp(h_t/2)');
box off; xlim([T_id(1)-1 T_id(end)+1]);

figure; 
hist(store_theta(:,1),50);
title('Histogram of posterior draws of psi');
box off;

% figure; 
% plot(psigrid,ppsi_mean); box off;
% title('Posterior density of psi');
% box off;

if cp_ml
    ml_BVAR_CSV_MA;
end



