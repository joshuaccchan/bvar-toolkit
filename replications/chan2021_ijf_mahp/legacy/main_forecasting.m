% This is the main run file for the forecasting exercise in Chan (2021)
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
%
% Sample period: 1959:Q1 - 2018:Q4
% Evaluation period: 1985:Q1 - end of sample
% Forecast horizons: 1- and 4-step-ahead

clear; clc;
nsim = 20000;
burnin = 100;
p = 4;

    % 1: Minnesota-type normal-gamma; 2: normal-gamma; 3: Minnesota;        
model = 2;

load 'macrodata_Q_2018Q4.csv';
data = macrodata_Q_2018Q4(1:238,:);
var_id = [1,2,22,23,35,37,57,58,59,76,81,83,95,120,138,144,145,147,148,152,160,161,245]; % n = 23
% var_id = [1,2,22,23,35,37,57,58,59,76,81,83,95,120,138,144,145,147,148,152,...
%    160,161,245,3,18,97,123,133,156,199]; % n = 30
Y0 = data(1:8,var_id);  % save the first 8 obs as the initial conditions
Y = data(9:end,var_id);
[T,n] = size(Y);
T0 = 91; % time index for 1974Q1

    % priors 
ah = zeros(n,1); Vh = 10*ones(n,1);
nuh0 = 5*ones(n,1); Sh0 = .01*ones(n,1).*(nuh0-1);
c01 = [1, 1/.04];      % prior for kappa1
c02 = [1, 1/.04^2];    % prior for kappa2
c03 = [1, 1/1];        % prior for kappa3
c04 = [1, 1/100];      % prior for kappa4
lam0_nu_psi = 1;       % prior for nu_psi:

yhat1 = zeros(T-T0,3*n+1);    % 1-quarter-ahead; [observed y, point forecast, log prelike, joint log prelike]
yhat4 = zeros(T-4-T0+1,3*n+1);  % 4-quarter-ahead           

switch model    
    case 1      
        model_name = 'BVAR-MNG';
        store_kappa_hat = zeros(T-T0,2);
        store_kappaCI = zeros(T-T0,2,2);       
    case 2
        model_name = 'BVAR-NG';
        store_kappa_hat = zeros(T-T0,1);
        store_kappaCI = zeros(T-T0,2);        
    case 3        
        model_name = 'BVAR-Minn';
        store_kappa_hat = zeros(T-T0,2);
        store_kappaCI = zeros(T-T0,2,2);       
end
disp(['Starting the recursive forecasting exercise for ' model_name '....']);

start_time = clock; 
for t = T0:T-1
    disp([ num2str(T-t+1) ' more loops to go... ' ] );     
    Yt = Y(1:t,:);    
    Tt = size(Yt,1);    
    tmpYt = [Y0(end-p+1:end,:); Yt];
    Zt = zeros(Tt,n*p); 
    for ii=1:p
        Zt(:,(ii-1)*n+1:ii*n) = tmpYt(p-ii+1:end-ii,:);
    end
    Zt = [ones(Tt,1) Zt];    
  
    if model == 1
        forecast_BVAR_MNG;
    elseif model == 2
        forecast_BVAR_NG;
    elseif model == 3
        forecast_BVAR_Minn;  
    end

        % store optimal hyperparameter values
    switch model    
        case 1
            store_kappa_hat(t-T0+1,:) = kappa_hat;    
            store_kappaCI(t-T0+1,:,:) = kappaCI;
        case 2
            store_kappa_hat(t-T0+1,:) = kappa_hat;    
            store_kappaCI(t-T0+1,:) = kappaCI;             
        case 3
            store_kappa_hat(t-T0+1,:) = kappa_hat;    
            store_kappaCI(t-T0+1,:,:) = kappaCI;           
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
fprintf(['RMSFEs for ' model_name ':'])
fprintf('\n'); 
fprintf('                      | 1-step-ahead, 4-step-ahead\n'); 
fprintf('GDP                   | %.1f, %.1f\n', RMSFE_1(1), RMSFE_4(1)); 
fprintf('Industrial production | %.1f, %.1f\n', RMSFE_1(3), RMSFE_4(3)); 
fprintf('Unemployment rate     | %.1f, %.1f\n', RMSFE_1(9), RMSFE_4(9)); 
fprintf('PCE inflation         | %.1f, %.1f\n', RMSFE_1(13), RMSFE_4(13)); 
fprintf('CPI inflation         | %.1f, %.1f\n', RMSFE_1(14), RMSFE_4(14)); 
fprintf('Fed funds rate        | %.1f, %.1f\n', RMSFE_1(16), RMSFE_4(16)); 

fprintf('\n'); 
fprintf(['Average of log predictive likelihoods for ' model_name ':'])
fprintf('\n'); 
fprintf('                      | 1-step-ahead, 4-step-ahead\n'); 
fprintf('GDP                   | %.1f, %.1f\n', aveprelike_1(1), aveprelike_4(1)); 
fprintf('Industrial production | %.1f, %.1f\n', aveprelike_1(3), aveprelike_4(3)); 
fprintf('Unemployment rate     | %.1f, %.1f\n', aveprelike_1(9), aveprelike_4(9)); 
fprintf('PCE inflation         | %.1f, %.1f\n', aveprelike_1(13), aveprelike_4(13)); 
fprintf('CPI inflation         | %.1f, %.1f\n', aveprelike_1(14), aveprelike_4(14)); 
fprintf('Fed funds rate        | %.1f, %.1f\n', aveprelike_1(16), aveprelike_4(16)); 


