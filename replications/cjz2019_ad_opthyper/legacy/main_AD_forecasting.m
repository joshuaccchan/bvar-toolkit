% This is the main run file for producing the forecasting results in
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
addpath(genpath('..'));
global  K_kn

nsims = 5000;
burnin = 0;

% 1: BVAR; 2: BVAR-3kappa; 3: BVAR-5kappa
model = 1;

load 'macrodata_Q_2018Q4_2.csv';
data = macrodata_Q_2018Q4_2;
var_id = [1,95,144,59,22,57,133,160,2,18,34,35,76,81,97,120,152,245]; % n = 18
Y0 = data(1:4,var_id);  % save the first 4 obs as the initial conditions
Y = data(5:end,var_id);
n = size(Y,2);
p = 4; % lag length
k = n*p+1;    

T0 = 95; % time index for 1984Q1
T = 234; % time index for 2018Q4

yhat1 = zeros(T-T0,3*n+1);    % 1-quarter-ahead; [observed y, point forecast, log prelike, joint log prelike]
yhat4 = zeros(T-4-T0+1,3*n+1);  % 4-quarter-ahead    
        
    % default hyperparameter values 
kappa_0 = [.05, 1, 100, 1, 1];


    % define a few things
options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
K_kn = commutation_matrix(k,n);
kappa_con_opt = kappa_0(1:3);
kappa_opt = kappa_0;
 
switch model    
    case 1      
        model_name = 'BVAR';
        store_kappa = zeros(T-T0,3);
    case 2
        model_name = 'BVAR-3kappa';
        store_kappa = zeros(T-T0,3);
    case 3
        model_name = 'BVAR-5kappa';
        store_kappa = zeros(T-T0,5);
end
disp(['Starting the recursive forecasting exercise for ' model_name '....']);

start_time = clock; 
for t = T0:T-1
    disp([ num2str(T-t+1) ' more loops to go... ' ] );
    
    Yt = Y(1:t,:);    
    Tt = size(Yt,1);    
    forecast_BVAR_NCP; 

        % store optimal hyperparameter values
    switch model    
      case 2
        store_kappa(t-T0+1,:) = kappa_con_opt;    
      case 3
        store_kappa(t-T0+1,:) = kappa_opt;    
    end    
   
        % store forecasts
    tmpmax = max(tmpyhat1(:,n+1:end));
    yhat1(t-T0+1,:)  = [Y(t+1,:) mean(tmpyhat1(:,1:n)) ...
        log(mean(exp(tmpyhat1(:,n+1:end)-tmpmax)))+tmpmax]; 
    if t<=T-4
        tmpmax = max(tmpyhat4(:,n+1:end));
        yhat4(t-T0+1,:)  = [Y(t+4,:) mean(tmpyhat4(:,1:n)) ...
            log(mean(exp(tmpyhat4(:,n+1:end)-tmpmax)))+tmpmax];     
    end
end

disp( ['Forecasting exercise takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
disp(' ' );
 
RMSFE_1 = sqrt(mean((yhat1(4:end,1:n) - yhat1(4:end,n+1:2*n)).^2))';
RMSFE_4 = sqrt(mean((yhat4(1:end,1:n) - yhat4(1:end,n+1:2*n)).^2))';
RMSFE = [RMSFE_1 RMSFE_4];
aveprelike_1 = mean(yhat1(4:end,2*n+1:end))';
aveprelike_4 = mean(yhat4(1:end,2*n+1:end))';
aveprelike = [aveprelike_1 aveprelike_4];

clc;
var_id = [1,8,11,15];
fprintf(['RMSFEs for ' model_name ' (in 1e3):'])
fprintf('\n'); 
fprintf('              | 1-step-ahead, 4-step-ahead\n');
fprintf('Real GDP      | %.4f, %.4f\n', RMSFE_1(var_id(1)), RMSFE_4(var_id(1))); 
fprintf('PCE inflation | %.4f, %.4f\n', RMSFE_1(var_id(3)), RMSFE_4(var_id(3))); 
fprintf('Fed funds rate| %.4f, %.4f\n', RMSFE_1(var_id(4)), RMSFE_4(var_id(4))); 
fprintf('Unemployment  | %.4f, %.4f\n', RMSFE_1(var_id(2)), RMSFE_4(var_id(2))); 
fprintf('\n'); 

fprintf(['Sum of log predictive likelihoods for ' model_name ':'])
fprintf('\n'); 
fprintf('              | 1-step-ahead, 4-step-ahead\n');
fprintf('All variables | %.1f, %.1f\n', aveprelike_1(end), aveprelike_4(end)); 
fprintf('Real GDP      | %.1f, %.1f\n', aveprelike_1(var_id(1)), aveprelike_4(var_id(1))); 
fprintf('PCE inflation | %.1f, %.1f\n', aveprelike_1(var_id(3)), aveprelike_4(var_id(3))); 
fprintf('Fed funds rate| %.1f, %.1f\n', aveprelike_1(var_id(4)), aveprelike_4(var_id(4))); 
fprintf('Unemployment  | %.1f, %.1f\n', aveprelike_1(var_id(2)), aveprelike_4(var_id(2))); 


