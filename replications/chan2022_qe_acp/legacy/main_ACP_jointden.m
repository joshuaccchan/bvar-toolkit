% This script reproduces the contour plot of the joint posterior density 
% of kappa1_tilde and kappa2_tilde in Chan (2022)
%
% This code is free to use for academic purposes only, provided that the 
% paper is cited as:
%
% Chan, J.C.C. (2022). Asymmetric Conjugate Priors for Large Bayesian VARs,
% Quantitative Economics, 13(3): 1145-1169
%
% This code comes without technical support of any kind. It is expected to
% reproduce the results reported in the paper. Under no circumstances will
% the authors be held responsible for any use (or misuse) of this code in
% any way.

clear; clc;
p = 5;          % if p > 8, need to change Y0 and Y below
addpath('./utility');
    % load data
data = xlsread('database_2019Q4.xlsx');
idx_ns = [1,2,4,5,10,11,12,13,15]; % index for variables in levels
Y0 = data(1:8,1:15);  % save the first 8 obs as the initial conditions
Y = data(9:end,1:15);
[T,n] = size(Y);
tmpY = [Y0(end-p+1:end,:); Y];
Z = zeros(T,n*p); 
for ii=1:p
    Z(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
end
Z = [ones(T,1) Z];
sig2 = get_resid_var(Y0,Y);
kappa = [.04,.04^2,1,100];

    % plot ml concour
[Kappa1,Kappa2] = meshgrid(0.01:.001:.2,.001:.0002:0.012);
store_lml = zeros(size(Kappa1,1),size(Kappa2,2));
for ii = 1:size(Kappa1,1)
    for ij = 1:size(Kappa2,2)
         store_lml(ii,ij) = ml_VAR_ACP(p,Y,Z,...
             prior_ACP_redu(n,p,[Kappa1(ii,ij),Kappa2(ii,ij),kappa(3),kappa(4)],sig2,idx_ns));
    end    
end
store_ml = exp(store_lml-max(max(store_lml)));
[ml_Sym,kappa_Sym] = get_OptSymKappa(Y0,Y,Z,p,'redu',idx_ns);
figure; subplot(1,2,1);
hold on
    plot(0:.01:.1,0:.01:.1,'k--');
    scatter(kappa_Sym(1),kappa_Sym(2),'filled');    
    contour(Kappa1,Kappa2,store_ml); 
hold off
box off; xlim([0 .1]); ylim([0 .1]);
xlabel('$\tilde{\kappa}_1$','interpreter','latex');
ylabel('$\tilde{\kappa}_2$','interpreter','latex');
pos=[0.2411,0.7095,-0.0804,-0.5286];
annotation('textarrow','Position',pos,'String','Symmetric prior');
subplot(1,2,2); 
hold on
    contour(Kappa1,Kappa2,store_ml); 
    scatter(.04,.0016,'s','filled');
hold off
xlim([0 .1]); ylim([0 .01]);
xlabel('$\tilde{\kappa}_1$','interpreter','latex');
ylabel('$\tilde{\kappa}_2$','interpreter','latex');
pos = [0.7536,0.1834,-0.03214,0.04286];
annotation('textarrow','Position',pos,'String','Subjective prior');
box off; colormap jet;
set(gcf,'Position',[200 200 800 400])
