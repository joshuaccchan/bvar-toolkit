clc; clear all;
nsim = 30000;
burnin = 5000; 
% WARNING: this might take a long time.
% For ease of reference, we have stored the four output .mat files in the 'results_mat' folder.
model=2; rev_option=0;
[Mean_kappa,Mean_Impact,Median_kappa,Median_Impact,Sig_mean,Sig_median] = func_main_SVAR_v2(model,rev_option,nsim,burnin);
save('fullsampleCS1.mat','Mean_kappa','Mean_Impact','Median_kappa','Median_Impact','Sig_mean','Sig_median');

model=2; rev_option=1;
[Mean_kappa,Mean_Impact,Median_kappa,Median_Impact,Sig_mean,Sig_median] = func_main_SVAR_v2(model,rev_option,nsim,burnin);
save('fullsampleCS2.mat','Mean_kappa','Mean_Impact','Median_kappa','Median_Impact','Sig_mean','Sig_median');

model=1; rev_option=0;
[Mean_kappa,Mean_Impact,Median_kappa,Median_Impact,Sig_mean,Sig_median] = func_main_SVAR_v2(model,rev_option,nsim,burnin);
save('fullsampleOI1.mat','Mean_kappa','Mean_Impact','Median_kappa','Median_Impact','Sig_mean','Sig_median');

model=1; rev_option=1;
[Mean_kappa,Mean_Impact,Median_kappa,Median_Impact,Sig_mean,Sig_median] = func_main_SVAR_v2(model,rev_option,nsim,burnin);
save('fullsampleOI2.mat','Mean_kappa','Mean_Impact','Median_kappa','Median_Impact','Sig_mean','Sig_median');