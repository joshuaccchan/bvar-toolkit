% bvar.forecast.realtime_loaddata - assembles the vintage-t real-time estimation panel.
%
%   [data_t,data_tpk] = bvar.forecast.realtime_loaddata(rt_data,nonrev_data, ...
%                           t,T0,tcode,var_type)
%
% rt_data     - struct of real-time vintage matrices, fields var1...varn, one
%               column per vintage
% nonrev_data - struct of non-revised series, same field names (read only for
%               var_type 4)
% t           - vintage index; the estimation panel ends at t-1
% T0          - period of the first vintage (the vintage column is t-T0+1 for
%               quarterly vintages, (t-T0)*3+3 for monthly ones)
% tcode       - n x 1 transformation code: 5 = 400*log-difference, 1 = level
% var_type    - n x 1 vintage/observation frequency code: 1 = quarterly
%               vintage, quarterly observations; 2 = quarterly vintage,
%               monthly observations averaged to quarters; 3 = monthly
%               vintage, monthly observations averaged to quarters;
%               4 = non-revised, monthly observations averaged to quarters
% data_t      - (t-2) x n estimation panel as seen in vintage t
% data_tpk    - 5 x n actual outturns from the LAST vintage, rows t-1:t+3
%
% See:
% Chan, J.C.C. (2020). Large Bayesian Vector Autoregressions. In: P. Fuleky (Eds),
% Macroeconomic Forecasting in the Era of Big Data, 95-125, Springer, Cham.
function [data_t,data_tpk] = realtime_loaddata(rt_data,nonrev_data,t,T0,tcode,var_type)
n = length(var_type);
data_t = zeros(t-2,n);
data_tpk = zeros(5,n); % actual observations from the last vintage
for i_var = 1:n
    var_name = ['var' num2str(i_var)];
    if var_type(i_var) == 1 % Q vin, Q obs;
        tmp_var_allvin = rt_data.(var_name);
        vin = t-T0+1; % vintage index
        tmp_var = tmp_var_allvin(:,vin);
        last_vin = tmp_var_allvin(:,end);
    elseif var_type(i_var) == 2 % Q vin, M obs
        tmp_var_allvin = rt_data.(var_name);
        vin = t-T0+1;
        tmp_var = tmp_var_allvin(:,vin);        
        n_q = floor(length(tmp_var)/3);
        tmp_var = mean(reshape(tmp_var(1:3*n_q),3,n_q))';
        last_vin = tmp_var_allvin(:,end);
        last_vin = mean(reshape(last_vin(1:3*n_q),3,n_q))';        
    elseif var_type(i_var) == 3 % M vin, M obs
        tmp_var_allvin = rt_data.(var_name);
        vin = (t-T0)*3+3;
        tmp_var = tmp_var_allvin(:,vin);        
        n_q = floor(length(tmp_var)/3);
        tmp_var = mean(reshape(tmp_var(1:3*n_q),3,n_q))';
        last_vin = tmp_var_allvin(:,end);
        last_vin = mean(reshape(last_vin(1:3*n_q),3,n_q))'; 
    elseif var_type(i_var) == 4 % non-revised, M obs        
        tmp_var = nonrev_data.(var_name);        
        n_q = floor(length(tmp_var)/3);
        tmp_var = mean(reshape(tmp_var(1:3*n_q),3,n_q))';
        last_vin = tmp_var;
    end    
    if tcode(i_var) == 5
        data_t(:,i_var) = 400 * log(tmp_var(2:t-1)./tmp_var(1:t-2));
        y_last_vin = 400 * log(last_vin(2:end)./last_vin(1:end-1));
        data_tpk(:,i_var) = y_last_vin(t-1:t+3); 
    elseif tcode(i_var) == 1
        data_t(:,i_var) = tmp_var(1:t-2);        
        data_tpk(:,i_var) = last_vin(t-1:t+3); 
    end
end

end