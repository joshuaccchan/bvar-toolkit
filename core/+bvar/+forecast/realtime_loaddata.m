% bvar.forecast.realtime_loaddata - assembles the vintage-t real-time estimation panel.
% Extracted 2026-09-01 (step 3): BYTE-IDENTICAL in chan2020_springer_largebvar/legacy/loaddata.m
% (canonical, this copy) and chan2020_jbes_kronecker/legacy/realtime_forecasts/loaddata.m.
% Function renamed loaddata -> realtime_loaddata. No unit test yet (needs the vintage structs);
% covered by the forecasting regression tests when drivers are functionized.
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