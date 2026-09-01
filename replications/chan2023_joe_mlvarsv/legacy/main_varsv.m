% This is the main run file for replicating the application in Chan (2022).
% In particular, this code computes the marginal likelihoods of 5 VARs:
% a standard homoskedastic VAR (VAR), a VAR with the common SV (VAR-CSV), 
% a VAR with Cholesky SV (VAR-SV), a VAR with factor SV (VAR-FSV), and an
% extension of VAR-SV with an outlier component (VAR-SVO)
%
% This code is free to use for academic purposes only, provided that the 
% paper is cited as:
%
% Chan, J.C.C. (2022). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, forthcoming.
%
% This code comes without technical support of any kind.  It is expected to
% reproduce the results reported in the paper. Under no circumstances will
% the authors be held responsible for any use (or misuse) of this code in
% any way.

addpath('./utility')  
clear; clc;
    % 1: VAR; 2: VAR-CSV; 3: VAR-SV; 4: VAR-FSV; 5: VAR-SVO
model = 2;
is_kappafixed = false;  % if true, kappa1 and kappa2 are fixed
is_kappasym = false;    % if true, kappa1 = kappa2
p = 4;          % if p > 8, need to change Y0 and Y below
r = 2;          % # of factors
nsim = 10000;   % # of MCMC iterations
burnin = 1000;  
M = 10000;      % # of importance sampling draws
cp_ml = 1;      % 1: compute marginal likelihood 
flag_marg = 2;  % 2: integrate out A/(A,Sig) and sig2
 
    % load data
data_all = load('macrodata_Q_2019Q4.csv'); % 1959Q1 to 2019Q4
% varid = [1,22,59,120,133,144,148]; % n = 7
varid = [1,22,59,120,133,144,148,2,35,57,81,95,152,160,245]; % n = 15
% varid = [1,22,59,120,133,144,148,2,35,57,81,95,152,160,245,...
%     3,18,23,37,58,76,83,97,123,138,145,147,157,161,199]; % n = 30

data = data_all(:,varid); 
Y0 = data(1:8,:);  % save the first 9 obs as the initial conditions
Y = data(9:end,:);

[T,n] = size(Y);    
y = reshape(Y',n*T,1);
tmpY = [Y0(end-p+1:end,:); Y];
X = zeros(T,n*p);
for i=1:p
    X(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end
X = [ones(T,1) X];    
k_beta = n*(n-1)/2;       % dimension of B0
k_alp = n^2*p + n;        % dimension of A
k = k_alp/n;

    % priors
kappa = .2^2; 
kappa1 = .2^2; kappa3 = 100; kappa4 = .2^2;
if is_kappasym
    kappa2 = kappa1; 
else
    kappa2 = (.2^2)^2;
end

switch model 
    case 1
        model_name = 'VAR-NCP';
        [Hyper.A0,Hyper.VA,Hyper.nu0,Hyper.S0] = prior_NCP(p,kappa,kappa3,Y0,Y);  
    case 2
        if is_kappafixed
            model_name = 'VAR-CSV-fixed-kappa';
        else
            model_name = 'VAR-CSV';
        end
        [Hyper.A0,Hyper.VA,Hyper.nu0,Hyper.S0] = prior_NCP(p,kappa,kappa3,Y0,Y);
        Hyper.c0 = [1,1/.2^2]; 
        Hyper.nuh = 3; Hyper.Sh = .1*(Hyper.nuh-1);
        Hyper.phi0 = .98; Hyper.Vphi = .05^2;        
    case 3        
        if is_kappafixed
            model_name = 'VAR-SV-fixed-kappa'; 
        elseif is_kappasym
            model_name = 'VAR-SV-sym';
        else
            model_name = 'VAR-SV'; 
        end
        [Hyper.alp0,Hyper.Valp,sig2_hat,U_hat] = prior_Minn(p,kappa1,kappa2,kappa3,Y0,Y);
        [Hyper.beta0,Hyper.Vbeta] = prior_B0(Y0,Y,kappa4);
        Hyper.c0 = [1,1/.2^2; 1,1/.2^2; 1,1]; 
        Hyper.nuh = 3*ones(n,1); Hyper.Sh = .1*(Hyper.nuh-1);
        Hyper.mu0 = zeros(n,1); Hyper.Vmu = 100*ones(n,1);
        Hyper.phi0 = .98*ones(n,1); Hyper.Vphi = .05^2*ones(n,1);        
    case 4 
        if is_kappafixed
            model_name = 'VAR-FSV-fixed-kappa'; 
        elseif is_kappasym
            model_name = 'VAR-FSV-sym';
        else
            model_name = 'VAR-FSV'; 
        end        
        [Hyper.alp0,Hyper.Valp,sig2_hat] = prior_Minn(p,kappa1,kappa2,kappa3,Y0,Y);
        Hyper.c0 = [1,1/.2^2; 1,1/.2^2;];
        Hyper.nuh = 3*ones(n+r,1); Hyper.Sh = .1*(Hyper.nuh-1);
        Hyper.mu0 = zeros(n+r,1); Hyper.Vmu = 100*ones(n+r,1);
        Hyper.phi0 = .98*ones(n+r,1); Hyper.Vphi = .05^2*ones(n+r,1);         
        Hyper.l0 = 0;
        Hyper.Vl = 1;     
    case 5
        if is_kappafixed
            model_name = 'VAR-SVO-fixed-kappa'; 
        elseif is_kappasym
            model_name = 'VAR-SVO-sym';
        else
            model_name = 'VAR-SVO'; 
        end
        [Hyper.alp0,Hyper.Valp,sig2_hat,U_hat] = prior_Minn(p,kappa1,kappa2,kappa3,Y0,Y);
        [Hyper.beta0,Hyper.Vbeta] = prior_B0(Y0,Y,kappa4);
        Hyper.c0 = [1,1/.2^2; 1,1/.2^2; 1,1]; 
        Hyper.nuh = 3*ones(n,1); Hyper.Sh = .1*(Hyper.nuh-1);
        Hyper.mu0 = zeros(n,1); Hyper.Vmu = 100*ones(n,1);
        Hyper.phi0 = .98*ones(n,1); Hyper.Vphi = .05^2*ones(n,1);    
        Hyper.p0a = 10/4; Hyper.p0b = (1-1/16)*40;
end       

% run MCMC + compute ML if cp_ml == true
disp(['Starting MCMC for ' model_name '...']);
switch model 
    case 1
        VAR_NCP;        
    case 2
        VAR_CSV;        
    case 3
        VAR_ARSV_redu; 
    case 4
        VAR_FSV;
    case 5
        VAR_ARSVO_redu; 
end

if cp_ml ~= 0
    fprintf(['log marginal likelihood of ' model_name ': %.1f \n'], lml); 
end

