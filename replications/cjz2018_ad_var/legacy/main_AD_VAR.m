% This is the main run file for conducting the prior sensitivity analysis
% in Chan, Jacobi and Zhu (2018).
%
% This code is free to use for academic purposes only, provided that the 
% paper is cited as:
%
% Chan, J.C.C., Jacobi, L. and Zhu, D. (2018). How Sensitive Are VAR 
% Forecasts to Prior Hyperparameters? An Automated Sensitivity Analysis,
% CAMA Working Paper 25/2018.
%
% This code comes without technical support of any kind.  It is expected to
% reproduce the results reported in the paper. Under no circumstances will
% the authors be held responsible for any use (or misuse) of this code in
% any way.

clear all; clc;
is_fullsample = 1; % 1: full sample; 0: short sample from 1989:Q3

Burn=500; Sample=5000; 

global T n  Y X YTY XTY XTX  nu_1 S_0 invB invVb K b_0 B_0 nu_0
folder = fileparts(which(mfilename)); 
addpath(genpath(folder));
% Data preprocessing
Y = xlsread('USdata_2017Q4.xlsx');
if is_fullsample
    Y = [Y(29:end, 2), Y(29:end, 4), Y(29:end, 3)];  % start from 1954:Q3
else
    Y = [Y(169:end, 2), Y(169:end, 4), Y(169:end, 3)]; % start from 1989:Q3
end
[T,n] =size(Y);
lag=2;
K=1+n*lag; 
d=n*K;
X=ones(T,K);
for i=3:T
   X(i,:)=[1,Y(i-1,:),Y(i-2,:)];
end    


%% Minnesota prior setup
ar_order = 4;
kepa_1 = 0.2^2; 
kepa_2 = 100;
kepa_3 = 1;
[b_0, B_0, v_0, S_0] = Min_Prior(Y, lag, ar_order, kepa_1, kepa_2, kepa_3);
K_qm = commutation_matrix(K, n);
I_p = eye(K);
f1 = kron(kron(eye(n), K_qm), eye(K));
ind1=(0:n-1)*K+1;
ind2=[2,3,4,9,10,11,16,17,18];
ind3=[5,6,7,12,13,14,19,20,21];
%% precompute
invB.v=B_0.v\speye(d);
invB.d=sparse(d_minverse(invB.v,B_0.d));
invVb.v=invB.v*b_0;
invVb.d=sparse(kron(b_0',eye(d))*invB.d);
nu_0=5; 
nu_1=nu_0+T;
XTX=X'*X;XTY=X'*Y;YTY=Y'*Y;
B=XTX\XTY;BigX=sparse(kron(eye(n),X)); BigXTY=sparse(kron(eye(n), XTY));
%% Simulation
% rng(10);
t=1:n;

id2=(t-1)*n+t;
N=Burn+Sample; forecast_period=20;
Z=randn(d,N);
%intialize sigma
iSigma.v=eye(n);
iSigma.d=zeros(n*n,3);
Ym.v=zeros(Sample,forecast_period,n);
Ym.d=zeros(forecast_period,n,3);
bm=zeros(K,n);
for g=1:N

   Vg.v=(invB.v+kron(iSigma.v,XTX))\speye(d);
   Vg.d=d_minverse(Vg.v,invB.d+f1*kron(iSigma.d,XTX(:)));
   
   
   bg.v=invVb.v+ reshape(XTY*iSigma.v,d,1);
   bg.d=invVb.d+BigXTY*iSigma.d; 
   bg.d=d_prod(Vg.v,Vg.d,bg.v,bg.d);
   bg.v=Vg.v*bg.v;

   C=chol(Vg.v)';
   b.v=bg.v+C*Z(:,g);
   b.d=bg.d+kron(Z(:,g)',eye(d))*d_cholasky(C,Vg.d);
   
   
   b.v=reshape(b.v,K,n);
   error.v=(Y-X*b.v); 
   error.d=-BigX*b.d;
   
   delta.v=S_0.v+error.v'*error.v;
   delta.d=S_0.d+ d_prod(error.v',d_trans(error.v,error.d),error.v,error.d);
   
   delta.v=delta.v\speye(n);
   delta.d=d_minverse(delta.v,delta.d);
   
   iSigma.v=wishrnd(delta.v,nu_1);
   L.v=chol(delta.v)';
   L.d=d_cholasky(L.v,delta.d);
   A=(L.v\speye(n))*iSigma.v;
   A2=iSigma.v*(L.v'\speye(n));
   
   iSigma.d=kron(A',eye(n))*L.d+kron(eye(n),A2)*d_trans(L.v,L.d);
   
   Sigma.v=iSigma.v\speye(n);  
   Sigma.d=d_minverse(Sigma.v,iSigma.d);
   
   if g>Burn
   ym=forecastVAR(ind1, ind2, ind3, forecast_period,b,Sigma);
   Ym.v(g-Burn,:,:)=ym.v;%+repmat(data0_mean,forecast_period,1);
   Ym.d=Ym.d+ym.d;
   Sigm.v(g-Burn,:)=diag(Sigma.v)';
   Sigm.d(g-Burn,:,:)=Sigma.d(id2,:);
   id=forecast_period*(g-Burn-1)+1:forecast_period*(g-Burn);
   forecastmean.v(id,:)=ym.vm;
   forecastmean.d(id,:,:)=ym.dm;
   bm=bm+b.v;
   end    
end 
bm=bm./Sample;
%Y=Y+repmat(data0_mean,T,1);
means=squeeze(mean(Ym.v,1));
Q_upper=squeeze(quantile(Ym.v,0.84,1));
Q_lower=squeeze(quantile(Ym.v,0.16,1));
Ym.d=Ym.d./Sample;
%% Sensitivity Output

for i=1:n   
    d_Qupper(i,:,:)=d_normal(Sigm.v(:,i),squeeze(Sigm.d(:,i,:)),Q_upper(:,i),forecastmean.v(:,i),...
                             squeeze(forecastmean.d(:,i,:)),forecast_period, Sample);   
                         
    d_QLower(i,:,:)=d_normal(Sigm.v(:,i),squeeze(Sigm.d(:,i,:)),Q_lower(:,i),forecastmean.v(:,i),...
                             squeeze(forecastmean.d(:,i,:)),forecast_period, Sample);      
end    

%% Forecast plots

plotforecast;

% %% Data plot
% tid = linspace(1955,2017.75,T)';
% figure
% subplot(3,1,1); plot(tid,Y(:,1),'linewidth',1);
% xlim([tid(1)-.5 tid(end)+.5]); ylim([3 11]); box off; title('Unemployment'); 
% subplot(3,1,2); plot(tid,Y(:,2),'linewidth',1);
% xlim([tid(1)-.5 tid(end)+.5]); box off; title('Interest Rate');
% subplot(3,1,3); plot(tid,Y(:,3),'linewidth',1);
% xlim([tid(1)-.5 tid(end)+.5]); ylim([-11 12]);  box off; title('GDP Growth');











