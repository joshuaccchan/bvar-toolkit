% This script estimates the BVAR-t-CSV model
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

%% prior
S0 = eye(n); nu0 = n+3;
nuh0 = 5; Sh0 = .01*(nuh0-1); 
rho0 = .9; Vrho = .2^2;
nuub = 100; %% upperbound for nu
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
store_lam = zeros(nsims,T);
store_theta = zeros(nsims,3);

counth = 0; countrho = 0; countnu = 0;
nugrid = linspace(2,nuub,700)';
store_pnu = zeros(700,1);

%% initialize the chain
h = zeros(T,1);
nu = 5;
rho = .8;
sigh2 = .1;
Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
lam = 1./gamrnd(nu/2,2/nu,T,1);


%% MCMC starts here
randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
disp('Starting MCMC for BVAR-t-CSV.... ');
start_time = clock;

for isim = 1:nsims + burnin
  
        %% sample Sig and A    
    iOm = sparse(1:T,1:T,exp(-h)./lam);
    XiOm = X'*iOm;
    KA = sparse(1:k,1:k,1./VA0) + XiOm*X;
    Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XiOm*shortY);
    Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + shortY'*iOm*shortY ...
        - Ahat'*KA*Ahat;
    Shat = (Shat+Shat')/2; % adjust for rounding errors
    Sig = iwishrnd(Shat,nu0+T);    
    CSig = chol(Sig,'lower');
    A = Ahat + (chol(KA,'lower')'\randn(k,n))*CSig'; 
    
        %% sample h
    U = shortY - X*A;
    tmp = (U/CSig');
    s2 = sum(tmp.^2,2)./lam;
    [h flag] = sample_h(s2,rho,sigh2,h,n);
    counth = counth + flag;
    
        %% sample lam
    U = shortY - X*A;
    tmp = (U/CSig');
    s2 = sum(tmp.^2,2)./exp(h);
    lam = 1./gamrnd((n+nu)/2,2./(s2+nu));
    
        %% sample nu        
    [nu flag fnu] = sample_nu(lam,nu,nuub);
    countnu = countnu + flag;    
    
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
    
    
    if isim > burnin
        isave = isim - burnin; 
        store_A = store_A + A;
        store_Sig = store_Sig + Sig;
        store_h(isave,:) = h';
        store_lam(isave,:) = lam';
        store_theta(isave,:) = [nu rho sigh2];
        
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
theta_mean = mean(store_theta)';
h_mean = mean(store_h)';
pnu_mean = store_pnu/nsims;
CSV_std_mean = mean(exp(store_h/2))';

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

T_id = linspace(1959,2013.75,T)';
figure; 
plot(T_id, CSV_std_mean); box off;
title('Posterior mean of exp(h_t/2)');
box off; xlim([T_id(1)-1 T_id(end)+1]);

if cp_ml
    ml_BVAR_t_CSV;
end
