This zip file contains the Matlab code for replicating the results 
in Chan, Koop, and Yu (2021). Specifically, it contains two parts:

Part I: full sample results: 

main_SVAR_fullsample.m, which uses func_main_SVAR_v2.m;
func_main_SVAR_v2 contains the main estimation programs: SVARSV_MH.m (OI-SV) and CS_MH.m (CS-SV);
after that, Fig89_plot_varcov_first2last2.m and Fig10_plot_B0_fullsample.m plot Figures 8-10.

Part II: forecasting results:

submain_forecasting_CS1.m, ..._CS2.m, ..._OI1.m, ..._OI2.m,
which use the main forecasting programs (basically estimation but adding a forecasting module):  
forecast_SVARSV_MH.m (OI-SV) and forecast_CS_MH.m (CS-SV).
(Note: in using submain_forecasting_... code, you have to input a "t", which represents a certain estimation point 
within our empirical data; in our implementation, we used cluster to help us get the full results, i.e., for all
t=sample_start_point,...sample_end_point. The full results are stored in the 'results_mat' folder.)
After that, Table3_forecasting.m generates Table 3.

This code is free to use for academic purposes only, and we would appreciate it if you could
cite our paper below:

Chan, J. C. C., Koop, G. and Yu, X. (2021). Large Order-Invariant 
Bayesian VARs with Stochastic Volatility. Working Paper.

This code comes without technical support of any kind. It is expected to
reproduce the results reported in the paper. Under no circumstances will
the author be held responsible for any use (or misuse) of this code in any way.

For questions, please contact Xuewen Yu: yu656@purdue.edu