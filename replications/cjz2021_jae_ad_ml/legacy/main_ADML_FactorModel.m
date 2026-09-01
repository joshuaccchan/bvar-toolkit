% This is the main run file for estimating the factor model and the 
% associated marginal likelihood in Chan, Jacobi and Zhu (2021). 
%
% This code is free to use for academic purposes only, provided that the 
% paper is cited as:
%
% Chan, J.C.C., Jacobi, L. and Zhu, D. (2021). An Automated Prior Robustness
% Analysis in Bayesian Model Comparison, Journal of Applied Econometrics, 
% forthcoming.
%
% This code comes without technical support of any kind.  It is expected to
% reproduce the results reported in the paper. Under no circumstances will
% the authors be held responsible for any use (or misuse) of this code in
% any way.

clear; clc;
% Determine where your m-file's folder is.
folder = fileparts(which(mfilename)); 
% Add that folder plus all subfolders to the path.
addpath(genpath(folder));
rng("default");
nsim = 10000;
burnin = 1000;
R = 10000;  % simulation size for CE method
r = 3;      % no of factors;

data = readmatrix('exchange_rates_2020.xlsx','Sheet','23 series','Range','B2:X1250');
var_id = [1,3,4,5,7,8,11,12,17,21];
returns = 100*log(data(2:end,var_id)./data(1:end-1,var_id));
[T,n] = size(returns);
Y = returns - repmat(mean(returns),T,1);
y = reshape(Y',T*n,1);
na = n*r - r*(r+1)/2;

    % initialize for storage
store_a.v = zeros(nsim,na);store_a.d=zeros(na,4,nsim);
%store_f.v = zeros(nsim,T*r);store_f.d=zeros(T*r,4,nsim);
store_Sig.v = zeros(nsim,n);store_Sig.d=zeros(n,4,nsim);
store_Omega.v = zeros(nsim,r);store_Omega.d=zeros(r,4,nsim);
store_invSig.v = zeros(nsim,n);store_invSig.d=zeros(n,4,nsim);
store_invOmega.v = zeros(nsim,r);store_invOmega.d=zeros(r,4,nsim);
store_lpa.v=zeros(nsim,1);store_lpa.d=zeros(nsim,4);
    % prior
a0 = 0; Va = 1; % aij iid N(a0,Va)
nusig2 = 3; Ssig2 = 1*(nusig2-1)*ones(n,1);
nuomega2 = 3; Somega2 = 1*(nuomega2-1)*ones(r,1);
prior = @(ax,s,o) prior_AD(ax,s,o,a0,Va,nusig2,Ssig2,nuomega2,Somega2);
prior2=@(ax,s,o) prior_AD_ordinates(ax,s,o,a0,Va,nusig2,Ssig2,nuomega2,Somega2);
     
rand('state', sum(100*clock) ); randn('state', sum(200*clock) );
disp(['Starting MCMC for the ' num2str(r) '-factor model.... ']);
    % initialize the Markov chain
% rng("default");
[invSig,invOmega,A,Ind_a,Ind,...
    K_rt,K_rn,K_rr,I_r,I_rr,E_r,Ka,Ia,Iaa,Ea]=Initialise_parameters(Y,na,n,r,4);
 idx=  1+(0:na-1)*(na+1);
Get_A=@(F,invSig)Sample_A(F,invSig,Y,idx, Va,a0,na,Ind,Ind_a,K_rt,K_rn,...
    Ka,Ea,Ia,Iaa);
Get_F=@(invSig,invOmega, A)Sample_F(invSig,invOmega,A,K_rn,K_rr,E_r,I_r,I_rr,r,T, Y);
Get_Sig=@(A,Ft)Sample_Sig(A,Ft,Y,Ssig2,T,nusig2,n);
Get_Omega=@(F)Sample_Omega(F,Somega2,nuomega2,T,r);  
In1=sparse(n*(0:n-1)+(1:n),1:n,ones(n,1));
In2=sparse(r*(0:r-1)+(1:r),1:r,ones(r,1));

for isim = 1:nsim+burnin    
    % sample f
    invSig.d =sparse(In1 * invSig.d);
    invSig.v=sparse(1:n,1:n,invSig.v);
    
    invOmega.d=sparse(In2*invOmega.d);
    invOmega.v=sparse(1:r,1:r,invOmega.v);
    F=Get_F(invSig,invOmega,A);
    
    % Sample A
    [A,a,Ft,K_a,mu_a,L]=Get_A(F,invSig);
    
    % sample Sig
    [Sig,invSig]=Get_Sig(A,Ft);
    
    % Sample Omega
    [Omega,invOmega]=Get_Omega(F);
    
    if (mod(isim, 1000) == 0)
        disp([num2str(isim) ' loops... '])
    end
    
    if isim>burnin
        i = isim-burnin;        
        store_a.v(i,:) = a.v';
        store_a.d(:,:,i)=a.d;
        store_invSig.v(i,:) = invSig.v';
        store_invSig.d(:,:,i)=invSig.d;
        store_Sig.v(i,:) = Sig.v';
        store_Sig.d(:,:,i)=Sig.d;
        store_invOmega.v(i,:) = invOmega.v';
        store_invOmega.d(:,:,i)=invOmega.d;
        store_Omega.v(i,:) = Omega.v';
        store_Omega.d(:,:,i)=Omega.d;
        l= lmvnpdfAD2(a,mu_a,K_a,L,na);
        store_lpa.v(i,1) = l.v;
        store_lpa.d(i,:)=l.d;
    end
    
end
disp('MCMC done');

disp('Computing ML using Chib.... ');
lpa=logmean(store_lpa,nsim);
ML_Chibs = SF_Chib_AD(lpa,store_a,store_Sig,store_Omega,...
    Y,prior2, Get_F,Get_Omega,Get_Sig, nusig2,nuomega2);
disp( 'ML computation via Chib done' );

disp('Computing ML using CE.... ');
ML_CE = SF_CE_AD(store_a,store_invSig,store_invOmega,Y,prior,R,Ind,Ind_a);
disp( 'ML computation via CE done' );

disp(' ')
fprintf('                                    CE          Chib\n');
fprintf('log-ML values:                      %.1f   %.1f\n', ML_CE.v, ML_Chibs.v);
fprintf('Derivative of log-ML w.r.t. kappa4:  %.1f     %.1f\n', ML_CE.d(1), ML_Chibs.d(1));
fprintf('Derivative of log-ML w.r.t. kappa5: %.1f     %.1f\n', ML_CE.d(2), ML_Chibs.d(2));
fprintf('Derivative of log-ML w.r.t. kappa6: %.0f     %.0f\n', ML_CE.d(3), ML_Chibs.d(3));
fprintf('Derivative of log-ML w.r.t. kappa7: %.1f    %.1f\n', ML_CE.d(4), ML_Chibs.d(4));

