% This is the main run file for estimating the VAR-t and the associated 
% marginal likelihood in Chan, Jacobi and Zhu (2021). 
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
global b_0 B_0 v_0 S_0 BigX n X Y
 % Determine where your m-file's folder is.
 folder = fileparts(which(mfilename)); 
% % Add that folder plus all subfolders to the path.
 addpath(genpath(folder));

lag = 2;
nu = 5; 
data = xlsread('USdata_2019Q4.xlsx','B2:E289');
var_id = [1,3];
Y0 = data(1:40,var_id);  % use the first 40 obs as presample
Y = data(41:end,var_id);
[T, n] = size(Y);
K = n*lag+1;
d = n*K;
X = ones(T,1);
for ii=1:lag
    X = [X,[Y0(end-ii+1:end,:);Y(1:end-ii,:)]];
end
X_t = X';
BigY=Y'; BigY=BigY(:);
q = n*(n*lag+1);
Kn=commutation_matrix(T,n);
%% Minnesota prior setup
kappa_1 = 0.4^2; 
kappa_2 = 100; 
kappa_3 = 1; 
k=3;

[b_0, B_0, v_0, S_0] = Min_Prior(Y0, lag, 4, kappa_1, kappa_2, kappa_3);
beta.v=(X_t*X)\(X_t*Y);BigX = kron(speye(n),X); BigXTY = kron(speye(n),X_t*Y);
beta.d=zeros(q,3); 
error.v=(Y-X*beta.v);error.d=-BigX*beta.d;
newnu=0.5*(n+nu);
Sigma.inv=eye(n);Sigma.invd=zeros(n*n,3);
%% The Gibbs loops
Sample = 10000; 
Burn = 1000; 
Z=randn(q,Sample+Burn);
G=gamrnd(newnu,1,T,Sample+Burn);

randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
disp('Starting MCMC for VAR-t.... ');
for i=1:Sample+Burn
    %sample lambda
    L=Getlambda(G(:,i),error,Sigma,nu,k, T);
    %sample beta
    BigS=Mtimes(X_t,L);
    C=Mtimes(BigS,X);
    B_g.v=B_0.inv+kron(BigS.v*X,Sigma.inv); B_g.d=B_0.invd+d_kron(C.v,C.d,Sigma.inv,Sigma.invd);
    B_g=Minverse(B_g);
    b_g.v=kron(BigS.v,Sigma.inv)*BigY;
    b_g.d=kron(BigY',speye(q))*d_kron(BigS.v,BigS.d,Sigma.inv, Sigma.invd);
    [beta,b_g, error]=d_beta2(B_g,b_g,Z(:,i),q);

    %sample Sigma
   delta=Minverse(Maddition(S_0,Mtimes(Mtrans(error,Kn),Mtimes(L,error)))); 
   Sigma=d_Sigma(delta,n,v_0+T);   
  
     if i>Burn
       ind=i-Burn; 
       Bm.v(ind,:)=beta.v';Bm.d(:,:,ind)=beta.d;
       Sigm.v(:,:,ind)=Sigma.v;Sigm.d(:,:,ind)=Sigma.d;
       Vb.v(:,:,ind)=B_g.v;Vb.d(:,:,ind)=full(B_g.d);
       bm.v(ind,:)=b_g.v';bm.d(:,:,ind)=full(b_g.d);
     end
end
disp('MCMC done');

disp('Computing ML for VAR-t using Chib.... ');
ML_Chibs = VAR_t_Chibs_AD(Bm,Sigm,bm,Vb,T,nu, K);
disp( 'ML computation via Chib done' );

disp('Computing ML for VAR-t using CE.... ');
ML_CE = VAR_t_CE_AD(Bm,Sigm,T, nu,K);
disp( 'ML computation via CE done' );

disp(' ')
fprintf('                                    CE          Chib\n');
fprintf('log-ML values for VAR-t:          %.1f   %.1f\n', ML_CE.v, ML_Chibs.v);
fprintf('Derivative of log-ML w.r.t. kappa1: %.1f     %.1f\n', ML_CE.d(1), ML_Chibs.d(1));
fprintf('Derivative of log-ML w.r.t. kappa2: %.2f     %.2f\n', ML_CE.d(2), ML_Chibs.d(2));
fprintf('Derivative of log-ML w.r.t. kappa3:  %.2f      %.2f\n', ML_CE.d(3), ML_Chibs.d(3));
