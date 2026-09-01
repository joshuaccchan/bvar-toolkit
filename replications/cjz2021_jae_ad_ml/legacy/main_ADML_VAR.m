% This is the main run file for estimating the VAR and the associated 
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

clc; clear;
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));
global T n  Y X YTY XTY XTX  nu_1 S_0 invB invVb K b_0 B_0 nu_0 BigX
lag = 2;
data = xlsread('USdata_2019Q4.xlsx','B2:E289');
var_id = [1,3];
Y0 = data(1:40,var_id);  % use the first 40 obs as presample
Y = data(41:end,var_id);
[T, n] = size(Y);
K=1+n*lag; 
d=n*K;
X = ones(T,1);
for ii=1:lag
    X = [X,[Y0(end-ii+1:end,:);Y(1:end-ii,:)]];
end

%% Minnesota prior setup
ar_order = 4;
kappa_1 = .4^2;
kappa_2 = 10^2;
kappa_3 = 1;
[b_0, B_0, nu_0, S_0] = Min_Prior(Y0, lag, ar_order, kappa_1, kappa_2, kappa_3);
f1 = kron(kron(speye(n), commutation_matrix(K, n)), speye(K));
%% precompute
invB.v=B_0.v\speye(d);
invB.d=sparse(d_minverse(invB.v,B_0.d));
invVb.v=invB.v*b_0.v;
invVb.d=sparse(kron(b_0.v',eye(d))*invB.d); 
nu_1=nu_0+T;
XTX=X'*X;XTY=X'*Y;YTY=Y'*Y;
B=XTX\XTY;BigX=kron(speye(n),X); BigXTY=kron(speye(n), XTY);

%% Simulation
Burn = 1000; Sample = 10000;
N=Burn+Sample; 
Z=randn(d,N);
%intialize sigma
Sigma.inv=eye(n);
Sigma.invd=zeros(n*n,3);

randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
disp('Starting MCMC for VAR.... ');
for g=1:N
    Vg.v=(invB.v+kron(Sigma.inv,XTX))\speye(d);
    Vg.d=d_minverse(Vg.v,invB.d+f1*kron(Sigma.invd,XTX(:)));
    bg.v=invVb.v+ reshape(XTY*Sigma.inv,d,1);
    bg.d=invVb.d+BigXTY*Sigma.invd;
    [b,bg, error]=d_beta(Vg,bg,Z(:,g),d);
    delta.v=S_0.v+error.v'*error.v;
    delta.d=S_0.d+ d_prod(error.v',d_trans(error.v,error.d),error.v,error.d);
    
    delta.v=delta.v\speye(n);
    delta.d=d_minverse(delta.v,delta.d);
    
    Sigma=d_Sigma(delta,n,nu_1);
    
    if g>Burn
        S.v(:,:, g-Burn)=Sigma.v;
        S.d(:,:,g-Burn)=Sigma.d;
        beta.v(g-Burn,:)=b.v;beta.d(:,:,g-Burn)=b.d;
        Varg.v(:,:,g-Burn)=Vg.v;Varg.d(:,:,g-Burn)=Vg.d;
        Bg.v(:,g-Burn)=bg.v; Bg.d(:,:,g-Burn)=bg.d;
    end
end 
disp('MCMC done');

disp('Computing ML for VAR using Chib.... ');
ML_Chibs = VAR_Chib_AD(beta,S,Varg,Bg);
disp( 'ML computation via Chib done' );

disp('Computing ML for VAR using CE.... ');
ML_CE= VAR_CE_AD(beta,S);
disp( 'ML computation via CE done' );

disp(' ')
fprintf('                                    CE          Chib\n');
fprintf('log-ML values for VAR:              %.1f   %.1f\n', ML_CE.v, ML_Chibs.v);
fprintf('Derivative of log-ML w.r.t. kappa1: %.1f     %.1f\n', ML_CE.d(1), ML_Chibs.d(1));
fprintf('Derivative of log-ML w.r.t. kappa2: %.2f     %.2f\n', ML_CE.d(2), ML_Chibs.d(2));
fprintf('Derivative of log-ML w.r.t. kappa3:  %.2f      %.2f\n', ML_CE.d(3), ML_Chibs.d(3));

    
figure;
tid = linspace(1948,2019.75,T)';
plot(tid, Y); box off;
xlim([tid(1)-1 tid(end)+2]);
legend('inflation','real GDP growth');
