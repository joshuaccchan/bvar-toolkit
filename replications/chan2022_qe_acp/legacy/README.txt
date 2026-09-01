This zip file contains Matlab code for replicating the empirical
results in Chan (2022). The main file to estimate the 6- and 15-variable
VAR identified with sign restrictions is 

main_ACP_apps.m 

The main file to reproduce the contour plot of the joint posterior density 
of kappa1_tilde and kappa2_tilde is 

main_ACP_jointden.m

The dataset database_2019Q4.xlsx consists of 15 US quarterly variables. 
The sample period is from 1985:Q1 to 2019:Q4. These variables are constructed
from raw time-series obtained from various sources, including the FRED 
database at the Federal Reserve Bank of St. Louis and the Federal Reserve
Bank of Philadelphia. The complete list of these time-series and their
sources are given in Appendix A of the paper.

This code is free to use for academic purposes only, provided that the 
paper is cited as:

Chan, J.C.C. (2022). Asymmetric Conjugate Priors for Large Bayesian VARs,
Quantitative Economics, 13(3): 1145-1169

This code comes without technical support of any kind. It is expected to
reproduce the results reported in the paper. Under no circumstances will
the author be held responsible for any use (or misuse) of this code in any way.