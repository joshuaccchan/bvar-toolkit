% GOLDEN-RUN VARIANT PATCH (2026-09-01): identical to legacy/main_BVAR.m except
% the model selector below is set to 4 instead of the shipped default 8
% (BVAR-CSV-t-MA), to capture marginal-likelihood goldens for all 8 models.
% This is the main run file for estimating the large Bayesian VARs 
% in Chan (2020). It also computes the corresponding marginal likelihoods.
%
% This code is free to use for academic purposes only, provided that the 
% paper is cited as:
%
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.
%
% This code comes without technical support of any kind.  It is expected to
% reproduce the results reported in the paper. Under no circumstances will
% the authors be held responsible for any use (or misuse) of this code in
% any way.

clear; clc;
% 1: BVAR; 2: BVAR-t; 3: BVAR-CSV; 4: BVAR-MA; 5: BVAR-t-CSV; 6: BVAR-t-MA;
% 7: BVAR-CSV-MA; 8: BVAR-CSV-t-MA
model = 4;
cp_ml = 1;      % 1: compute marginal likelihood 
p = 4;          % if p > 4, need to change Y0 and shortY below
nsims = 30000;
burnin = 5000;

%% load data
load 'data_Q.csv'; % 1959Q1 to 2013Q4
data = data_Q(:,[1:3 6:15 17 19:24]); % The list of variables is given in Appendix B
Y0 = data(1:4,:);  % save the first 4 obs as the initial conditions
shortY = data(5:end,:);
[T,n] = size(shortY);    
k = n*p+1;         % # of coefficients in each equation

switch model
    case 1 
        BVAR;
    case 2        
        BVAR_t;
    case 3
        BVAR_CSV;      
    case 4
         BVAR_MA;
    case 5
        BVAR_t_CSV;
    case 6
        BVAR_t_MA;
    case 7
        BVAR_CSV_MA;
    case 8
        BVAR_CSV_t_MA;
end    