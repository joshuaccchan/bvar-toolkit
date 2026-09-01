% This is the main run file for producing the full sample estimation 
% results in Chan (2019).
%
% This code is free to use for academic purposes only, provided that the 
% paper is cited as:
%
% Chan, J.C.C. (2019). Asymmetric conjugate priors for large Bayesian VARs,
% CAMA Working Papers 51/2019
%
% This code comes without technical support of any kind. It is expected to
% reproduce the results reported in the paper. Under no circumstances will
% the authors be held responsible for any use (or misuse) of this code in
% any way.

clear; clc;
p = 4;          % # of lags; if p > 8, need to change Y0 and Y below
nsim = 10000;   % # of posterior draws
is_contour = false; % construct contour plot for kappa1 and kappa2; true or false

    % load data
load 'macrodata_Q_2018Q4.csv';
data = macrodata_Q_2018Q4(1:238,:);
var_id = [1,2,18,22,34,35,57,59,76,81,95,97,120,123,133,138,145,148,152,160,245]; % n = 21

Y0 = data(1:8,var_id);  % save the first 8 obs as the initial conditions
Y = data(9:end,var_id);
[T,n] = size(Y);
tmpY = [Y0(end-p+1:end,:); Y];
Z = zeros(T,n*p); 
for ii=1:p
    Z(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
end
Z = [ones(T,1) Z];
kappa = [.04, .04^2, 1, 100]; % subjective prior

    % obtain optimal parameters
[ml_opt,kappa_opt] = get_OptKappa(Y0,Y,Z,p,[.04,.04]);  % asymmetric prior
[ml_Sym,kappa_Sym] = get_OptSymKappa(Y0,Y,Z,p);  % symmetric prior

fprintf('         | Symmetric prior | Asymmetric prior:')
fprintf('\n'); 
fprintf('kappa_1  | %.3f           | %.3f \n', kappa_Sym(1), kappa_opt(1)); 
fprintf('kappa_2  | %.3f           | %.3f \n', kappa_Sym(2), kappa_opt(2)); 
fprintf('log-ML   | %.0f           | %.0f  \n', ml_Sym, ml_opt); 
fprintf('\n'); 

    % sample from the posterior
[store_alp,store_beta,store_Sig] = sample_ThetaSig(Y0,Y,p,kappa_opt,nsim);

    % convert draws to the reduced form
[store_Btilde,store_Sigtilde] = getReducedForm(store_alp,store_beta,store_Sig);
    % the VAR coefficients are stacked by row - see the paper
Btilde_hat = mean(store_Btilde)'; 
Sigtilde_hat = squeeze(mean(store_Sigtilde));

if is_contour
    % contour plot for kappa1 and kappa2
[Kappa1,Kappa2] = meshgrid(0.25:.01:.65,.002:.0002:0.02);
store_lml = zeros(size(Kappa1,1),size(Kappa2,2));
for ii = 1:size(Kappa1,1)
    for ij = 1:size(Kappa2,2)
         store_lml(ii,ij) = ml_VAR_ACP(p,Y,Y0,Z,[Kappa1(ii,ij),Kappa2(ii,ij),kappa(3),kappa(4)]);
    end    
end
store_ml = exp(store_lml-max(max(store_lml)));

figure; subplot(1,2,1);
hold on
    plot(0:.01:.6,0:.01:.6,'k--');
    scatter(.039,.039,'filled');
    scatter(.04,.0016,'s','filled');
    contour(Kappa1,Kappa2,store_ml); 
hold off
box off; xlim([0 .6]);
subplot(1,2,2); contour(Kappa1,Kappa2,store_ml); 
box off; colormap jet;
end