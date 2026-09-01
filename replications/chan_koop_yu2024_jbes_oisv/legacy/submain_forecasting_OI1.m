% This is the main run file for the forecasting exercise in Chan, Koop, and Yu (2021)
% using OI-1
% Sample period: 1959:03 - 2019:12
% Evaluation period: 1970:03 - end of sample
% Forecast horizons: 1-12 step-ahead

%clear; clc;
addpath('./utility')
nsim = 10000;
burnin = 500;
p = 13;

% 1: Minnesota-type Horseshoe SVARSV_redu; 2: Minnesota-type Horseshoe CCCM_redu;
model = 1;

data = load('FRED_MD_20vars.csv'); %4--IP, 6--unemployment, 12--PCE inflation, 13--FEDFUNDS % start from 1959/03
var_id1 = [4,6,12,13]; % n = 20
var_id2 = [1:3,5,7:11,14:20]; % first four variables: GDP growth, unemployment, PCE inflation, FEDFUNDS rate
var_id = [var_id1, var_id2]; %var_id = flip(var_id);
Y0 = data(1:24,var_id);  % save the first 24 obs as the initial conditions
Y = data(25:end,var_id);
[T,n] = size(Y);
T0 = 108; % time index for 1970:02

yhat_mean = zeros(3*n+3,12);    % 1~12-month-ahead; [observed y, point forecast, log prelike, joint log prelike]
yhat_median = zeros(3*n+3,12);    % 1~12-month-ahead; [observed y, point forecast, log prelike, joint log prelike]

switch model
    case 1
        model_name = 'fore-SVARSV-MH';
        store_kappa_hat = zeros(T-T0,2);
        store_kappaCI = zeros(T-T0,2,2);
        % priors
        Hyper.nuh = 3*ones(n,1); Hyper.Sh = .05*ones(n,1).*(Hyper.nuh-1);
        Hyper.phi0 = .95*ones(n,1); Hyper.Vphi = .05^2*ones(n,1);
        Hyper.B0 = eye(n); Hyper.VB0 = 1*ones(n); % prior variances of each elements
    case 2
        model_name = 'fore-CS-MH';
        store_kappa_hat = zeros(T-T0,2);
        store_kappaCI = zeros(T-T0,2,2);
        % priors
        Hyper.nuh = 3*ones(n,1); Hyper.Sh = .05*ones(n,1).*(Hyper.nuh-1);
        Hyper.phi0 = .95*ones(n,1); Hyper.Vphi = .05^2*ones(n,1);
        kappa = [.1,.1,1,100]; % kappa(1): own lag; kappa(2): other lag; kappa(4): intercepts; kappa(3): alpha
        k_beta = n^2*p + n; Hyper.beta0 = zeros(k_beta,1);
        Hyper.alp0 = zeros(n*(n-1)/2,1); Hyper.Valp = 1*ones(n*(n-1)/2,1);
        Hyper.mu0 = zeros(n,1); Hyper.Vmu = 100*ones(n,1);
end
disp(['Starting the recursive forecasting exercise for ' model_name '....']);

start_time = clock;
%for t = T0:T-1
disp([ num2str(T-t+1) ' more loops to go... ' ] );
Yt = Y(1:t,:);
Tt = size(Yt,1);
tmpYt = [Y0(end-p+1:end,:); Yt];
Xt = zeros(Tt,n*p);
for ii=1:p
    Xt(:,(ii-1)*n+1:ii*n) = tmpYt(p-ii+1:end-ii,:);
end
Xt = [ones(Tt,1) Xt];

if model == 1
    forecast_SVARSV_MH;
elseif model == 2
    forecast_CS_MH;
end

% store forecasts
for tt=1:12
    tmpmax = max(tmpyhat(:,n+1:end,tt));
    if t<=T-tt
        yhat_median(:,tt)  = [Y(t+tt,:) median(tmpyhat(:,1:n,tt)) ...
            log(mean(exp(tmpyhat(:,n+1:end,tt)-tmpmax)))+tmpmax];
    else
        yhat_median(:,tt)  = [NaN(1,n) median(tmpyhat(:,1:n,tt)) ...
            log(mean(exp(tmpyhat(:,n+1:end,tt)-tmpmax)))+tmpmax];
    end
end

for tt=1:12
    tmpmax = max(tmpyhat(:,n+1:end,tt));
    if t<=T-tt
        yhat_mean(:,tt)  = [Y(t+tt,:) mean(tmpyhat(:,1:n,tt)) ...
            log(mean(exp(tmpyhat(:,n+1:end,tt)-tmpmax)))+tmpmax];
    else
        yhat_mean(:,tt)  = [NaN(1,n) mean(tmpyhat(:,1:n,tt)) ...
            log(mean(exp(tmpyhat(:,n+1:end,tt)-tmpmax)))+tmpmax];
    end
end

disp( ['Forecasting exercise takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
disp(' ' );

%save(['.\results\' model_name],'yhat','store_kappa','store_kappaCI');




