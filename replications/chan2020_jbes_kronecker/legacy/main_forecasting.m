% This is the main run file for producing the forecasting results in 
% Chan (2020).
%
% This code is free to use for academic purposes only, provided that the 
% paper is cited as:
%
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.
%
% This code comes without technical support of any kind. It is expected to
% reproduce the results reported in the paper. Under no circumstances will
% the authors be held responsible for any use (or misuse) of this code in
% any way.
% 
% sample period: 1964Q1 - 2015Q4
% evaluation period: 1975Q1 - 2013Q4
% forecast horizons: 0, 1, 2, 4

addpath('./realtime_forecasts');
clear; clc;
nsims = 2000;
burnin = 100;

% 1: BVAR; 2: BVAR-t; 3: BVAR-CSV; 4: BVAR-MA; 5:BVAR-t-CSV; 6:BVAR-t-MA; 
% 7: BVAR-CSV-MA; 8: BVAR-t-CSV-MA; 9: BVAR-SV-small
model = 3;

n = 20; % number of variables
p = 4; % lag length
k = n*p+1;
T0 = 41; % time index for 1974Q1 
T = 208; % time index for 20154Q4
var_small = [1 8 12 14]; % real GDP growth, unemployment, PCE inflation, Fed Funds rate

% load all vintages
rt_data.var1 = xlsread('ROUTPUTQvQd.xlsx','AI70:HA283');
rt_data.var2 = xlsread('RCONQvQd.xlsx','AI70:HA283');
rt_data.var3 = xlsread('rinvbfQvQd.xlsx','AI70:HA283');
rt_data.var4 = xlsread('rinvresidQvQd.xlsx','AI70:HA283');
rt_data.var5 = xlsread('RNXQvQd.xlsx','AI71:HA283');
rt_data.var6 = xlsread('npiQvQd.xlsx','AI70:HA283');
rt_data.var7 = xlsread('iptMvMd.xlsx','EF542:YI1184');    % monthly vintage, monthly obs
rt_data.var8 = xlsread('rucQvMd.xlsx','AI209:HA848');     % quarterly vintage, monthly obs
rt_data.var9 = xlsread('employMvMd.xlsx','DG302:XJ944');  % monthly vintage, monthly obs
rt_data.var10 = xlsread('hMvMd.xlsx','AD206:UG848');      % monthly vintage, monthly obs
rt_data.var11 = xlsread('hstartsMvMd.xlsx','BW206:VX848');% monthly vintage, monthly obs
rt_data.var12 = xlsread('pconQvQd.xlsx','AI70:HA283');  
rt_data.var13 = xlsread('pimpQvQd.xlsx','AI70:HA283');  
nonrev_data.var14 = xlsread('FEDFUNDS.xls','B129:B770');  % monthly obs
nonrev_data.var15 = xlsread('GS1.xls','B144:B785');       % monthly obs
nonrev_data.var16 = xlsread('GS10.xls','B144:B785');      % monthly obs
nonrev_data.var17 = xlsread('BAAFFM.xls','B129:B770');    % monthly obs
nonrev_data.var18 = xlsread('ISM-MAN_PMI.xls','B196:B837'); % monthly obs
nonrev_data.var19 = xlsread('ISM-MAN_NEWORDERS.xls','B196:B837'); % monthly obs
nonrev_data.var20 = xlsread('SP500.xlsx','B62:B705');    % monthly obs
tcode = [5 5 5 5 1 5 5 1 5 5 5 5 5 1 1 1 1 1 1 5]';
    % 1: Q vin, Q obs; 2:Q vin, M obs; 3: M vin, M obs; 4: non-revised, M obs
var_type = [1 1 1 1 1 1 3 2 3 3 3 1 1 4 4 4 4 4 4 4]; 

yhat0 = zeros(T-T0+1,3*4+1);   % nowcast; [observed y (store only 4 variables), point forecasts, log prelike, joint log prelike] 
yhat1 = zeros(T-1-T0+1,3*4+1); % 1-step-ahead
yhat2 = zeros(T-2-T0+1,3*4+1); % 2-step-ahead
yhat4 = zeros(T-4-T0+1,3*4+1); % 4-step-ahead
   
%% prior
if model <=8
    S0 = eye(n); nu0 = n+3;
    psi0 = 0; Vpsi = 1;
    lpri_psi = @(x) -.5*(x-psi0)^2/Vpsi -1e10*(x<-.99 || x>.99);
    nuub = 100;       %% upperbound for nu
    nuh0 = 5; Sh0 = .01*(nuh0-1); 
    rho0 = .9; Vrho = .2^2;
    kappa1 = .2^2; kappa2 = 100;
    A0 = zeros(k,n);
    VA0 = zeros(k,1);
    sig2 = zeros(n,1);
end

switch model    
    case 1      
        model_name = 'BVAR';
    case 2
        model_name = 'BVAR-t';
    case 3
        model_name = 'BVAR-CSV';
    case 4
        model_name = 'BVAR-MA';
    case 5
        model_name = 'BVAR-t-CSV';
    case 6
        model_name = 'BVAR-t-MA';
    case 7
        model_name = 'BVAR-CSV-MA';
    case 8
        model_name = 'BVAR-t-CSV-MA';
    case 9
        model_name = 'BVAR-SV-small';
end
disp(['Starting the recursive forecasting exercise for ' model_name '....']);

start_time = clock; 
for t = T0:T
    disp([ num2str(T-t+1) ' more periods to go... ' ] );
    
    [data_t,data_tpk] = loaddata(rt_data,nonrev_data,t,T0,tcode,var_type);
    
        % check if the lastest data are missing 
    is_last_miss = (sum(isnan(data_t(end,:)),2)>0);
    if is_last_miss 
        data_t = data_t(1:end-1,:);
    end
        % discard past history with missing values
    ind = sum(isnan(data_t),2);
    ind_last = find(ind > 0, 1, 'last' );
    if ~isempty(ind_last)
        data_t = data_t(ind_last+1:end,:);
    end
    Y0 = data_t(1:4,:);
    shortYt = data_t(5:end,:);   
    Tt = size(shortYt,1);
    
    if model<= 8
        % construct VA0
    tmpY = [Y0(end-p+1:end,:); shortYt];
    for i=1:n
        Z = [ones(Tt,1) tmpY(4:end-1,i) tmpY(3:end-2,i) tmpY(2:end-3,i)...
            tmpY(1:end-4,i)];
        tmpb = (Z'*Z)\(Z'*tmpY(5:end,i));
        sig2(i) = mean((tmpY(5:end,i)-Z*tmpb).^2);
    end
    for i=1:k
        l = ceil((i-1)/n); 
        idx = mod(i-1,n); % variable index
        if idx==0
            idx = n;
        end        
        if i==1 % intercept
            VA0(1) = kappa2;     
        else
            VA0(i) = kappa1/(l^2*sig2(idx));        
        end   
    end
    X = zeros(Tt,n*p); 
    for i=1:p
        X(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
    end
    X = [ones(Tt,1) X];    
    end
  
        %% run the forecast model here
    switch model    
        case 1      
            forecast_BVAR;
        case 2
            forecast_BVAR_t;
        case 3
            forecast_BVAR_CSV;        
        case 4
            forecast_BVAR_MA;
        case 5
            forecast_BVAR_t_CSV;
        case 6
            forecast_BVAR_t_MA;
        case 7
            forecast_BVAR_CSV_MA;
        case 8       
            forecast_BVAR_t_CSV_MA;
        case 9
            forecast_BVAR_SV_small;
    end    
    tmpmax = max(tmpyhat0(:,[n+var_small,2*n+1]));
    yhat0(t-T0+1,:)  = [data_tpk(1,var_small) mean(tmpyhat0(:,var_small)) ...
        log(mean(exp(tmpyhat0(:,[n+var_small,2*n+1])-tmpmax)))+tmpmax]; 
    if t<=T-1
        tmpmax = max(tmpyhat1(:,[n+var_small,2*n+1]));
        yhat1(t-T0+1,:)  = [data_tpk(2,var_small) mean(tmpyhat1(:,var_small)) ...
            log(mean(exp(tmpyhat1(:,[n+var_small,2*n+1])-tmpmax)))+tmpmax];     
    end
    if  t<=T-2
        tmpmax = max(tmpyhat2(:,[n+var_small,2*n+1]));
        yhat2(t-T0+1,:) = [data_tpk(3,var_small) mean(tmpyhat2(:,var_small)) ...
            log(mean(exp(tmpyhat2(:,[n+var_small,2*n+1])-tmpmax)))+tmpmax];   
    end
    if  t<=T-4
        tmpmax = max(tmpyhat4(:,[n+var_small,2*n+1]));
        yhat4(t-T0+1,:) = [data_tpk(5,var_small) mean(tmpyhat4(:,var_small)) ...
            log(mean(exp(tmpyhat4(:,[n+var_small,2*n+1])-tmpmax)))+tmpmax];   
    end
end

disp( ['Forecasting exercise takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
disp(' ' );

RMSFE_0 = sqrt(mean((yhat0(5:end,1:4) - yhat0(5:end,5:8)).^2))'; % evaluation period: 1975Q1-2015Q4
RMSFE_1 = sqrt(mean((yhat1(4:end,1:4) - yhat1(4:end,5:8)).^2))';
RMSFE_2 = sqrt(mean((yhat2(3:end,1:4) - yhat2(3:end,5:8)).^2))';
RMSFE_4 = sqrt(mean((yhat4(:,1:4) - yhat4(:,5:8)).^2))';
RMSFE = [RMSFE_0 RMSFE_1 RMSFE_2 RMSFE_4];
aveprelike_0 = mean(yhat0(5:end,9:end))';
aveprelike_1 = mean(yhat1(4:end,9:end))';
aveprelike_2 = mean(yhat2(3:end,9:end))';
aveprelike_4 = mean(yhat4(:,9:end))';
aveprelike = [aveprelike_0 aveprelike_1 aveprelike_2 aveprelike_4];

save(model_name,'yhat0','yhat1','yhat2','yhat4');

clc;
fprintf(['Average of log predictive likelihoods for ' model_name ':\n'])
fprintf('             | real GDP, unemployment, PCE, Fed funds rate, joint\n'); 
fprintf('nowcast      | %.2f, %.2f, %.2f, %.2f, %.1f\n', aveprelike_0); 
fprintf('1-step-ahead | %.2f, %.2f, %.2f, %.2f, %.1f\n', aveprelike_1); 
fprintf('2-step-ahead | %.2f, %.2f, %.2f, %.2f, %.1f\n', aveprelike_2); 
fprintf('4-step-ahead | %.2f, %.2f, %.2f, %.2f, %.1f\n', aveprelike_4); 
fprintf('\n'); 

fprintf(['RMSFE for ' model_name ':\n'])
fprintf('             | real GDP, unemployment, PCE, Fed funds rate\n'); 
fprintf('nowcast      | %.2f, %.2f, %.2f, %.2f\n', RMSFE_0); 
fprintf('1-step-ahead | %.2f, %.2f, %.2f, %.2f\n', RMSFE_1); 
fprintf('2-step-ahead | %.2f, %.2f, %.2f, %.2f\n', RMSFE_2); 
fprintf('4-step-ahead | %.2f, %.2f, %.2f, %.2f\n', RMSFE_4); 


