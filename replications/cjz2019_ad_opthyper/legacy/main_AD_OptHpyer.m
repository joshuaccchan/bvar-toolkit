% This is the main run file for producing the full sample results in
% Chan, Jacobi and Zhu (2019)
%
% This code is free to use for academic purposes only, provided that the 
% paper is cited as:
%
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019.
%
% This code comes without technical support of any kind. It is expected to
% reproduce the results reported in the paper. Under no circumstances will
% the authors be held responsible for any use (or misuse) of this code in
% any way.

clear; clc;
addthisfiles
global  K_kn
p = 4;          % if p > 4, need to change Y0 and Y below

load 'macrodata_Q_2018Q4_2.csv';
data = macrodata_Q_2018Q4_2(1:238,:);
var_id = [1,95,144,59,22,57,133,160,2,18,34,35,76,81,97,120,152,245]; % n = 18
Y0 = data(1:4,var_id);  % save the first 4 obs as the initial conditions
Y = data(5:end,var_id);
[T,n] = size(Y);
k = n*p+1;    
kappa_0 = [.05, 1, 100, 1, 1];   % default hyperparameter values for kappa1,...,kappa4
options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
K_kn = commutation_matrix(k,n);

    % evaluate the lml at the default hyperparameter values
[A0,VA,nu0,S0] = prior_NCP(p,kappa_0,Y0,Y);
iVA = VA\speye(k);
tmpY = [Y0(end-p+1:end,:); Y];
Z = zeros(T,n*p); 
for i=1:p
    Z(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end
Z = [ones(T,1) Z];
KA = iVA + Z'*Z;
A_hat = KA\(iVA*A0 + Z'*Y);
S_hat = S0 + A0'*iVA*A0 + Y'*Y - A_hat'*KA*A_hat;
S_hat = (S_hat+S_hat')/2;
ml_0 = ml_VAR_NCP(VA,S0,nu0,KA,S_hat,T);

    % optimize kappa1 - kappa5
kappa_inl = kappa_0;
[kappa_opt,ml_opt,~,output] = fmincon(@(kappa)FiveKappa(kappa,Y,Y0,T,n,k,p),...
    kappa_inl,[],[],[],[],zeros(5,1),[Inf;10;Inf;Inf;Inf],[],options);
ml_opt = - ml_opt;

    % optimize kappa1 - kappa3
[kappa_con_opt,ml_con_opt,~,output_con] = fmincon(@(kappa)ThreeKappa(kappa,kappa_0(1,4:5),Y,Y0,T,n,k,p),...
    kappa_inl(1,1:3),[],[],[],[],zeros(3,1),[Inf;10;Inf],[],options);
ml_con_opt = -ml_con_opt;

fprintf('\n'); 
fprintf('                  | baseline | optimize kappa1-kappa3 | optimize kappa1-kappa5:\n'); 
fprintf('kappa1            | %.3f    | %.3f                  | %.3f\n', kappa_0(1), kappa_con_opt(1), kappa_opt(1)); 
fprintf('kappa2            | %.1f      | %.1f                    | %.1f\n', kappa_0(2), kappa_con_opt(2), kappa_opt(2)); 
fprintf('kappa3            | %.1f    | %.1f                   | %.1f\n', kappa_0(3), kappa_con_opt(3), kappa_opt(3)); 
fprintf('kappa4            | %.1f      | %.1f                    | %.1f\n', kappa_0(4), kappa_0(4), kappa_opt(4)); 
fprintf('kappa5            | %.1f      | %.1f                    | %.1f\n', kappa_0(5), kappa_0(5), kappa_opt(5)); 
fprintf('log-ML            | %.0f    | %.0f                  | %.0f\n', ml_0, ml_con_opt, ml_opt); 



