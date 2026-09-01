% This is the main run file for estimating the hybrid TVP-VAR in Chan (2022)
%
% This code is free to use for academic purposes only, provided that the 
% paper is cited as:
%
% Chan, J.C.C. (2022). Large Hybrid Time-Varying Parameter VARs, Journal of
% Business and Economic Statistics, forthcoming
%
% This code comes without technical support of any kind. It is expected to
% reproduce the results reported in the paper. Under no circumstances will
% the authors be held responsible for any use (or misuse) of this code in
% any way.

clear; clc;
addpath('./utility');
nsim = 50000; 
burnin = 1000;

% load data
p = 2;
data = xlsread('macrodata_Q_2018Q4.csv','A1:IN238'); % 1959Q1-2018Q4
% var_id = [1,95,59]; % n = 3
var_id = [1,95,59,144,22,133];   % n = 6
% var_id = [1,95,59,144,22,133,160,2,18,23,35,57,76,97,120,123,135,138,148,161]; % n = 20      
Y0 = data(1:4,var_id);  % save the first 4 obs as the initial conditions
Y = data(5:end,var_id);

[T,n] = size(Y);
y = reshape(Y',T*n,1);
k_alp = n*(n-1)/2;       % dimension of the impact matrix
k_beta = n^2*p + n;      % number of VAR coefficients

is_gamfixed = false;   % whether gam is fixed or estimated
is_kappafixed = false; % whether kappa is fixed or estimated

    % prior
kappa = [.4 .04^2];
sig2 = get_resid_var_v2(Y0,Y);
[C,idx_kappa1,idx_kappa2] = get_C(n,p,sig2);
[Valp,Vbeta] = getVtheta(idx_kappa1,idx_kappa2,kappa,C,sig2);
nuh0 = 3*ones(n,1); Sh0 = .1*(nuh0-1).*ones(n,1);    % priors for Sigh
Vsigbeta = .01^2*ones(k_beta,1);     % priors for Sigbeta
Vsigbeta(1:n*p+1:end) = .1^2;        % intercepts
Vsigalp = .01^2*ones(k_alp,1);       % priors for Sigalp
c01 = [1, 1/.04];      % prior for kappa1:
c02 = [1, 1/.04^2];    % prior for kappa2:
ah = zeros(n,1); Vh = 10*ones(n,1);  % priors for h0
ap = .5*ones(1,2); bp = .5*ones(1,2);

    % compute and define a few things
tmpY = [Y0(end-p+1:end,:); Y];
X2 = zeros(T,n*p); 
for ii=1:p
    X2(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
end
X2 = [ones(T,1) X2];

    % initialize for storage
store_alp = zeros(nsim,T,k_alp);
store_beta = zeros(nsim,T,k_beta);
store_h = zeros(nsim,T,n);
store_Sigbeta = zeros(nsim,k_beta);
store_Sigalp = zeros(nsim,k_alp);
store_Sigh = zeros(nsim,n);
store_beta0 = zeros(nsim,k_beta);
store_alp0 = zeros(nsim,k_alp);
store_h0 = zeros(nsim,n);
store_p0 = zeros(nsim,n,2);
store_gam = zeros(nsim,n*2);
store_kappa = zeros(nsim,2);
store_lpostgam = zeros(nsim,4);

    % initialize the Markov chain
if is_gamfixed
    gam = [0,0; [zeros(n-1,1),ones(n-1,1)]];
    % gam = zeros(n,2);
else
    gam = zeros(n,2);
end
U = zeros(T,n);
lp_gam = zeros(n,4);
initialize;

    % MCMC starts here
randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
disp('Starting MCMC for the hybrid TVP-VAR...');
start_time = clock;
 
for isim = 1:nsim + burnin  
        % sample states and parameters equation by equation     
    [Valp,Vbeta] = getVtheta(idx_kappa1,idx_kappa2,kappa,C,sig2);     
    for ii=1:n            
        ki = n*p+1+ii-1;
        idx_a0 = (ii-1)*(ii-2)/2+1; idx_a1 = ii*(ii-1)/2; % index for alp
        idx_b0 = (ii-1)*k_beta/n+1; idx_b1 = ii*k_beta/n; % index for beta
        Xi = [X2 -Y(:,1:ii-1)];
        Yi = Y(:,ii);
        hi = h(:,ii);
        gami = gam(ii,:);
        p0i = p0(ii,:)';        
        SigBeta = reshape(Sigbeta,k_beta/n,n);
        Sigbetai = SigBeta(:,ii);
        Sigalpi = Sigalp(idx_a0:idx_a1);
        
            % sample gami and thetai
        Sigthetai = [Sigbetai; Sigalpi];
        thetai0 = [beta0(idx_b0:idx_b1); alp0(idx_a0:idx_a1)];        
        [gami,thetai,Ui,tilde_thetai,lp_gami] = sample_gam_thetai_ver2(Xi,...
            thetai0,Sigthetai,Yi,hi,p0i,n*p+1,is_gamfixed,gami);        
        gam(ii,:) = gami;
        lp_gam(ii,:) = lp_gami;
        Thetai = reshape(thetai,ki,T)';
        beta(:,idx_b0:idx_b1) = Thetai(:,1:k_beta/n);
        alp(:,idx_a0:idx_a1) = Thetai(:,k_beta/n+1:end);
        U(:,ii) = Ui;
          
            % sample hi
        hi = sample_SVRW(log(Ui.^2),hi,Sigh(ii),h0(ii));
        h(:,ii) = hi;
        
            % sample thetai0 and Sigthetai
        tilde_Thetai = reshape(tilde_thetai,ki,T)';
        Vmui = [Vbeta(idx_b0:idx_b1); Valp(idx_a0:idx_a1); Vsigbeta(idx_b0:idx_b1); Vsigalp(idx_a0:idx_a1)];        
        Wi = [Xi Xi.*tilde_Thetai];
        WiSig  = Wi'*sparse(1:T,1:T,exp(-hi));
        Kmui = sparse(1:2*ki,1:2*ki,1./Vmui) + WiSig*Wi;
        mui_hat = Kmui\(WiSig*Yi);
        mui = mui_hat + chol(Kmui,'lower')'\randn(2*ki,1);      
        beta0(idx_b0:idx_b1) = mui(1:k_beta/n);
        alp0(idx_a0:idx_a1) = mui(k_beta/n+1:ki);
        Sigbeta(idx_b0:idx_b1) = mui(ki+1:ki+k_beta/n).^2;
        Sigalp(idx_a0:idx_a1) = mui(k_beta/n+ki+1:end).^2;        
    end    
   
        % sample h0
    Kh0 = sparse(1:n,1:n,1./Sigh + 1./Vh);
    h0_hat = Kh0\(ah./Vh + h(1,:)'./Sigh);
    h0 = h0_hat + chol(Kh0,'lower')'\randn(n,1);
    
        % sample Sigh
    eh = h - [h0';h(1:T-1,:)];
    Sigh = 1./gamrnd(nuh0+T/2,1./(Sh0 + sum(eh.^2)'/2));    
       

        % sample kappa1 and kappa2
    if ~is_kappafixed
        tmpc1 = sum(beta0(idx_kappa1).^2./C(idx_kappa1));
        tmpc2 = sum(beta0(idx_kappa2).^2./C(idx_kappa2));
        kappa(1) = gigrnd(c01(1)-n*p/2,2*c01(2),tmpc1,1);
        kappa(2) = gigrnd(c02(1)-(n-1)*n*p/2,2*c02(2),tmpc2,1);            
    end
    
        % sample p0
    if ~is_gamfixed        
        p0 = betarnd(ap+gam,bp+1-gam);
    end  
    
        % evaluate posterior of gamma
    lpostgam(1) = sum(lp_gam(:,1));                   % HYB-(0,0)
    lpostgam(2) = lp_gam(1,1) + sum(lp_gam(2:end,2)); % HYB-(0,1)
    lpostgam(3) = sum(lp_gam(:,3));                   % HYB-(1,0)
    lpostgam(4) = lp_gam(1,3)+sum(lp_gam(2:end,4));   % HYB-(1,1)
     
    if isim > burnin
        isave = isim - burnin;
        store_h(isave,:,:) = h;
        store_alp(isave,:,:) = alp;
        store_beta(isave,:,:) = beta;
        store_Sigbeta(isave,:) = Sigbeta';
        store_Sigalp(isave,:) = Sigalp';
        store_Sigh(isave,:) = Sigh';
        store_beta0(isave,:) = beta0';
        store_alp0(isave,:) = alp0';
        store_h0(isave,:) = h0';
        store_p0(isave,:,:) = p0;
        store_gam(isim-burnin,:) = reshape(gam',1,n*2);
        store_kappa(isave,:) = kappa';
        store_lpostgam(isave,:) = lpostgam;
    end
    
    if (mod(isim, 1000) == 0)
        disp([num2str(isim) ' loops... ']);
    end     
end

disp( ['MCMC takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
disp(' ' );

% compulte SDDR
maxtmp = max(store_lpostgam);
lpostgam = log(mean(exp(store_lpostgam-repmat(maxtmp,nsim,1)))) + maxtmp;
lprigam_c = n*sum(gammaln(ap+bp)-gammaln(ap+bp+1)-gammaln(ap)-gammaln(bp));
lprigam = zeros(1,4);
lprigam(1) = lprigam_c + n*sum(gammaln(ap)+gammaln(bp+1));
lprigam(2) = lprigam_c + n*(gammaln(ap(1)+1)+gammaln(bp(1))+gammaln(ap(2))+gammaln(bp(2)+1));
lprigam(3) = lprigam_c + n*(gammaln(ap(1))+gammaln(bp(1)+1)+gammaln(ap(2)+1)+gammaln(bp(2)));
lprigam(4) = lprigam_c + n*sum(gammaln(ap+1)+gammaln(bp));
lBF = lprigam - lpostgam;
gam_mode = reshape(get_gammode(store_gam),2,n)';
gam_hat = reshape(mean(store_gam),2,n)';

disp('Equation: Posterior mean of gamma_i^{\beta}, Posterior mean of gamma_i^{\alpha}');
eqtext = ["real GDP", "PCE inflation", "Unemployment", "Fed funds rate",...
	"Industrial production index", "Real average hourly earnings in manufacturing"...
    "M1", "Real PCE", "Real disposable personal income", "Industrial production: final products", ...
	"All employees: total nonfarm", "Civilian employment", ...
    "Nonfarm business section: hours of all persons", "GDP deflator",...
    "CPI", "PPI", "Nonfarm business sector: real compensation per hour",...
    "Nonfarm business section: real output per hour",...
    "10-year treasury constant maturity rate", "M2"];
for jj=1:n
    tmp = append(eqtext(jj), ': %.2f, %.2f\n');    
    fprintf(tmp, gam_hat(jj,1), gam_hat(jj,2));
end
disp(' ' );

disp('Log Bayes factors of the proposed hybrid TVP-VAR against: ');
fprintf('HYB(1,1): %.0f\n', lBF(4));
fprintf('HYB(1,0): %.0f\n', lBF(3));
fprintf('HYB(0,1): %.0f\n', lBF(2));
fprintf('HYB(0,0): %.0f\n', lBF(1));


