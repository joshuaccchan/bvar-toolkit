% This is the main run file for estimating the large Bayesian VAR with the
% Minnesota-type normal-gamma prior in Chan (2021)
%
% This code is free to use for academic purposes only, provided that the 
% paper is cited as:
%
% Chan, J.C.C. (2021). Minnesota-Type Adaptive Hierarchical Priors for 
% Large Bayesian VARs, International Journal of Forecasting, forthcoming
%
% This code comes without technical support of any kind.  It is expected to
% reproduce the results reported in the paper. Under no circumstances will
% the authors be held responsible for any use (or misuse) of this code in
% any way.

clear; clc;
    % 1: Minnesota-type normal-gamma; 2: normal-gamma; 3: Minnesota
model = 1; 
nsim = 10000; burnin = 1000;
p = 4;          % if p > 8, need to change Y0 and Y below
    % load data
load 'macrodata_Q_2018Q4.csv';
data = macrodata_Q_2018Q4(1:238,:);
var_id = [1,2,22,23,35,37,57,58,59,76,81,83,95,120,138,144,145,147,148,152,160,161,245]; % n = 23
% var_id = [1,2,22,23,35,37,57,58,59,76,81,83,95,120,138,144,145,147,148,152,...
%     160,161,245,3,18,97,123,133,156,199]; % n = 30
Y0 = data(1:8,var_id);  % save the first 8 obs as the initial conditions
Y = data(9:end,var_id);

[T,n] = size(Y);
tmpY = [Y0(end-p+1:end,:); Y];
Z = zeros(T,n*p); 
for ii=1:p
    Z(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
end
Z = [ones(T,1) Z];
k_beta = n^2*p+n;
k_alp = n*(n-1)/2;
sig2 = get_resid_var(Y0,Y);
[C,idx_kappa1,idx_kappa2] = get_C(n,p,sig2);

    % priors
ah = zeros(n,1); Vh = 10*ones(n,1);
nuh0 = 5*ones(n,1); Sh0 = .01*ones(n,1).*(nuh0-1);
c01 = [1, 1/.04];      % prior for kappa1/kappa
c02 = [1, 1/.04^2];    % prior for kappa2
lam0_nu_psi = 1;       % prior for nu_psi

switch model
    case 1 
        BVAR_MNG;
    case 2        
        BVAR_NG;
    case 3
        BVAR_Minn;   
end  