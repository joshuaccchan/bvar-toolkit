% GOLDEN-RUN VARIANT PATCH (2026-09-01): identical to legacy/main_forecasting.m
% except the model selector is 1 instead of the shipped default 8 (BVAR-CSV-t-MA),
% to capture recursive-forecasting goldens for all 8 models.
% This is the main run file for producing the forecasting results in 
% Chan (2020).
%
% This code is free to use for academic purposes only, provided that the 
% paper is cited as:
%
% Chan, J.C.C. (2020). Large Bayesian Vector Autoregressions. In: P. Fuleky (Eds),
% Macroeconomic Forecasting in the Era of Big Data, 95-125, Springer, Cham
%
% This code comes without technical support of any kind. It is expected to
% reproduce the results reported in the paper. Under no circumstances will
% the authors be held responsible for any use (or misuse) of this code in
% any way.
% 
% sample period: 1964Q1 - 2015Q4
% evaluation period: 1975Q1 - 2015Q4
% forecast horizons: nowcasts (0), 1-quarter-ahead (1) 
clear; clc;
nsims = 5000;
burnin = 100;
% 1: small BVAR; 2:BVAR + Minnesota prior; 3: BVAR + natural conjugate prior
% 4: BVAR + independent prior; 5: BVAR + SSVS prior; 6: BVAR-CSV; 
% 7: BVAR-CSV-t; 8: BVAR-CSV-t-MA
model = 1;

n = 20; % number of variables
p = 4; % lag length
k = n*p+1;
T0 = 41; % time index for 1974Q1 
T = 208; % time index for 20154Q4
var_core = [1 7 8 12]'; % real GDP growth, IP, unemployment, PCE rate

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
nonrev_data.var18 = xlsread('ISM-MAN_PMI.csv','B196:B837'); % monthly obs
nonrev_data.var19 = xlsread('ISM-MAN_NEWORDERS.xls','B196:B837'); % monthly obs
nonrev_data.var20 = xlsread('SP500.xlsx','B62:B705');    % monthly obs
tcode = [5 5 5 5 1 5 5 1 5 5 5 5 5 1 1 1 1 1 1 5]';
    % 1: Q vin, Q obs; 2:Q vin, M obs; 3: M vin, M obs; 4: non-revised, M obs
var_type = [1 1 1 1 1 1 3 2 3 3 3 1 1 4 4 4 4 4 4 4]; 

if model == 1
    yhat0 = zeros(T-T0+1,3*4+1);  
    yhat1 = zeros(T-1-T0+1,3*4+1);    
else
    yhat0 = zeros(T-T0+1,3*n+1);    % nowcast; [observed y, point forecast, log prelike, joint log prelike]
    yhat1 = zeros(T-1-T0+1,3*n+1);  % 1-quarter-ahead    
end
    
    
%% prior
psi0 = 0; Vpsi = 1;
lpri_psi = @(x) -.5*(x-psi0)^2/Vpsi -1e10*(x<-.99 || x>.99);
nuub = 50;       %% upperbound for nu
nuh0 = 5; Sh0 = .01*(nuh0-1); 
rho0 = .9; Vrho = .2^2;
q = .5;
if model == 1 || model == 2 || model == 4 || model == 5 
    c1 = .2^2; c2 = .1^2; c3 = 100;
elseif model == 3 || model == 6 || model == 7 || model == 8
    c1 = .2^2; c2 = 100;
end

switch model    
    case 1
        model_name = 'BVAR-small';        
    case 2
        model_name = 'BVAR-Minn';
    case 3
        model_name = 'BVAR-NCP';
    case 4
        model_name = 'BVAR-IP';
    case 5
        model_name = 'BVAR-SSVS';
    case 6
        model_name = 'BVAR-CSV';
    case 7
        model_name = 'BVAR-CSV-t';
     case 8
        model_name = 'BVAR-CSV-t-MA';
end
disp(['Starting the recursive forecasting exercise for ' model_name '....']);

start_time = clock; 
for t = T0:T
    disp([ num2str(T-t+1) ' more loops to go... ' ] );
    
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
      
        %% run the forecast model here
    switch model    
        case 1
            forecast_BVAR_small;
        case 2      
            forecast_BVAR_Minn;
        case 3
            forecast_BVAR_NCP;
        case 4
            forecast_BVAR_IP; 
        case 5
            forecast_BVAR_SSVS;
        case 6
            forecast_BVAR_CSV;
        case 7
            forecast_BVAR_CSV_t;
        case 8       
            forecast_BVAR_CSV_t_MA;
    end    
    if model == 1
        data_tp0 = data_tpk(1,var_small);
        data_tp1 = data_tpk(2,var_small);
    else
        data_tp0 = data_tpk(1,:);
        data_tp1 = data_tpk(2,:);
    end    
    tmpmax = max(tmpyhat0(:,n+1:end));
    yhat0(t-T0+1,:)  = [data_tp0 mean(tmpyhat0(:,1:n)) ...
        log(mean(exp(tmpyhat0(:,n+1:end)-tmpmax)))+tmpmax]; 
    if t<=T-1
        tmpmax = max(tmpyhat1(:,n+1:end));
        yhat1(t-T0+1,:)  = [data_tp1 mean(tmpyhat1(:,1:n)) ...
            log(mean(exp(tmpyhat1(:,n+1:end)-tmpmax)))+tmpmax];     
    end
end

disp( ['Forecasting exercise takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
disp(' ' );

if model==1        
    RMSFE_0 = sqrt(mean((yhat0(5:end,1:n) - yhat0(5:end,n+1:2*n)).^2))'; % evaluation period: 1975Q1-2015Q4
    RMSFE_1 = sqrt(mean((yhat1(4:end,1:n) - yhat1(4:end,n+1:2*n)).^2))';
    aveprelike_0 = mean(yhat0(5:end,2*n+1:end))';
    aveprelike_1 = mean(yhat1(4:end,2*n+1:end))';    
else
    RMSFE_0 = sqrt(mean((yhat0(5:end,var_core) - yhat0(5:end,n+var_core)).^2))'; % evaluation period: 1975Q1-2015Q4
    RMSFE_1 = sqrt(mean((yhat1(4:end,var_core) - yhat1(4:end,n+var_core)).^2))';
    aveprelike_0 = mean(yhat0(5:end,2*n+var_core))';
    aveprelike_1 = mean(yhat1(4:end,2*n+var_core))';    
end
RMSFE = [RMSFE_0 RMSFE_1];
aveprelike = [aveprelike_0 aveprelike_1];

save(['.\' model_name],'yhat0','yhat1');

clc;
fprintf(['RMSFE for ' model_name ':'])
fprintf('\n'); 
fprintf('                  | GDP, IP, urate, PCE\n'); 
fprintf('now-casts         | %.1f, %.1f, %.1f, %.1f\n', RMSFE(1,1), RMSFE(2,1), RMSFE(3,1), RMSFE(4,1)); 
fprintf('1-quarter-ahead   | %.1f, %.1f, %.1f, %.1f\n', RMSFE(1,2), RMSFE(2,2), RMSFE(3,2), RMSFE(4,2)); 
fprintf('\n'); 

fprintf(['ALPL for ' model_name ':'])
fprintf('\n'); 
fprintf('                  | GDP, IP, urate, PCE\n'); 
fprintf('now-casts         | %.1f, %.1f, %.1f, %.1f\n', aveprelike(1,1),...
    aveprelike(2,1), aveprelike(3,1), aveprelike(4,1)); 
fprintf('1-quarter-ahead   | %.1f, %.1f, %.1f, %.1f\n', aveprelike(1,2),...
    aveprelike(2,2), aveprelike(3,2), aveprelike(4,2)); 
fprintf('\n'); 

