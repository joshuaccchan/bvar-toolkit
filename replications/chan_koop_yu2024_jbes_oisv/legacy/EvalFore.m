% This function computes the Diebold-Mariano
% statistic, provide by Roque Montero.

function results = EvalFore(mfore,typets,yf)
%-------------------------------------------------------------------------%
% Matlab 9.0
% Autor: Roque Montero
% Date: 17/Dec/2016
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
% Description: Computes RMSE, ratio of RMSE, DM test, small sample correction 
% of the DM test (Harvey et all 1997), Giacomini and White (2006) test and
% the Pesaran and Timmermann (1994) test. 
% 1./ Note: Benchmark models is located in the first column.
% 2./ Note: ratio RMSE < 1 indicates alternative model outperform benchmark. 
% For tests, negative value indicates that alternative model outperforms 
% benchmark model.
% Inputs:
%   mfore1 : Forecast errors. Each column is a diff. model.
%   typets : (1) one sided test; (2) two sided test (Default: 1);
%   yf     : Evaluation window for Pesaran and Timmermann (1994) test.
%
% Outputs:
%   results : 
%   -.RMSE   : Root MSE for each model.
%   -.rMSE   : Ratio MSE with respect to the first model mfore(:,:,1).
%   -.DM     : Diebold and Mariano test and Pvalue.
%   -.HR     : Harvey, Leybourne, and Newbold test and Pvalue.
%   -.GW     : Giacomini and White test and Pvalie.
%   -.DA     : Directional accuracy test and Pvalue.
%   -.table  : Table with results.
%
% Index.
% 1. Test.
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%    
% 1. Test
if exist('typets','var') == 0
    typets = 1;
end
% Defining variables for computation of RMSE.
T = size(mfore,1);
M = size(mfore,2);

% MSE
efore2 = mfore.^(2);
MSE   = sum(efore2)' / T;
% Ratios of MSE.
rMSE = zeros(M-1,1);
for i0 = 2:M
    rMSE(i0-1) = MSE(i0)/ MSE(1);
end
results.RMSE = MSE.^(0.5);
results.rMSE = rMSE;

% Diebold and Mariano (1995) and Harvey, Leybourne, and Newbold 
% (HLN, 1997) Forecast Accuracy and Giacomini and White (2006).
maxlag = round(12*(T/100)^(0.25));
% Matrix to store results DM.
DMtest = zeros(M-1,1);
DMPval = zeros(M-1,1);
% Matrix to store results HR.
Smlpcorr = (T + 1 - 2*maxlag + maxlag*(1 - maxlag)/T)/T;
HR     = zeros(M-1,1);
HRPval = zeros(M-1,1);
% Matrix to store results GW.
GW     = zeros(M-1,1);
GWPval = zeros(M-1,1);

% Statistics and Pvalues.
id = 1;
for i0 = 2:M
    % Diebold and Mariano (1995)
    % Square diff.
    dif2 = efore2(:,i0) - efore2(:,1);
    % Mean of dif2 (dvar).
    DM = mean(dif2);
    % Variance: HAC estimator of variance var(dvar).
    dmvar = sqrt(NWhac(dif2,maxlag));
    % Test statistic
    DMtest(id) = sqrt(T)*(DM / dmvar);
    
    % Harvey, Leybourne, and Newbold (HLN, 1997).
    HR(id) = DMtest(id) .* sqrt(Smlpcorr);

    % Giacomini and White (2006).
    GW(id) = sqrt(T)*(mean(dif2) / sqrt(NWhac(dif2,maxlag)));
        
    % Pvalue DM, HR and GW.
    if typets == 2
        DMPval(id) = 2*(1-normcdf(abs(DMtest(id))));
        HRPval(id) = 2*(1-tcdf(abs(HR(id)),T-1));
    elseif typets == 1
        DMPval(id) = normcdf(DMtest(id));
        HRPval(id) = tcdf(HR(id),T-1);    
    end   
    GWPval(id) = 1-cdf('chi2',GW(id)^2,1);
    
    % Next model.
    id = id + 1;
    clear dif2 DM dmvar;
end

% Accuracy direction test Pesaran and Timmermann (1994)
% Statistics and Pvalues.
if nargin == 3
    % Matrix to store results.
    yfore= repmat(yf,1,M) - mfore;
    Iaux = yfore.*repmat(yf,1,M);
    % Sign true data.
    Sg_yf = double(gt(yf,0));
    % Sign of forecast.
    Sg_yff= double(gt(yfore,0)); 
    % Counting righ predictions.
    Signs = double(gt(Iaux,0));
    % Computing Statistics for sign predictions.
    % Effective data, forecast, dirrection forecast.
    Py = sum(Sg_yf);   
    Px = sum(Sg_yff,1)';
    Phat = sum(Signs,1)'; 
    Py = Py / T;
    Px = Px / T;
    % Succes ratio.
    SR  = Phat / T;
    % Succes ratio under independence(SRI).
    SRI = Py.*Px + (1 - Py).*(1 - Px);  
    % Compute variance for p and SR.
    vSR    = (T^(-1)).*(SRI.*(1 - SRI));
    vSRI   = ((T.^(-1))*(2*Py - 1)^(2)).*Px.*(1 - Px) ...
           + ((T^(-1)).*(2*Px - 1).^(2))*(Py*(1 - Py)) ...
           + (4*(T^(-2))*Py).*((1 - Py)*Px.*(1 - Px));
    % Computing the DA statistic.DA statistic.
    DA = (SR - SRI) ./ ((vSR - vSRI).^(1/2));
    Pvl= 1 - normcdf(DA,0,1);
    % Results.
    results.DA = [DA Pvl];
end

% Resluts.
results.DM = [DMtest DMPval];
results.HR = [HR HRPval];
results.GW = [GW GWPval];

% Print results
temp_1 = num2cell([results.RMSE(1) NaN(1,7)]);
temp_2 = num2cell((2:M)');
if nargin == 2
    labels = {'Model' 'RMSE' 'rMSE' 'DM' 'Pval' 'HR' 'Pval' 'GW' 'Pval'};
    temp_3 = num2cell([results.RMSE(2:end,1) rMSE results.DM results.HR results.GW]);
    Results = ([labels; {' '} temp_1; temp_2 temp_3]);
elseif nargin == 3
    labels = {'Model' 'RMSE' 'rMSE' 'DM' 'Pval' 'HR' 'Pval' 'GW' 'Pval' 'DA' 'Pval'};
    temp_3 = num2cell([results.RMSE(2:end,1) rMSE results.DM results.HR results.GW results.DA(2:end,:)]);
    Results = ([labels; {' '} temp_1 num2cell(results.DA(1,:)); temp_2 temp_3]);
end
results.table = Results;
%display(Results);
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%
%-------------------------------------------------------------------------%

end

%%% --- Below are some functions built in the main function EvalFore.m ---
function Sigma = NWhac(Z,nlag)
T = length(Z(:,1));
Sigma = Z'*Z / T;
for i = 1:nlag
    Z1 = Z((i+1):end, :);
    Z2 = Z(1:(end-i), :);
    Gamma_lag = Z1' * Z2 / T;
    Sigma = Sigma + (1- i/(nlag+1)) * (Gamma_lag + Gamma_lag');
end
end