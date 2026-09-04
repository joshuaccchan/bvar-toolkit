% bvar.sv.init_approx1N - crude 1-component log-chi2 approximation used to initialize SV paths.
% Extracted 2026-09-01 (step 3, zero-risk core): verified identical (modulo comments/whitespace) across chan2023_joe_mlvarsv/legacy/utility/getARh_approx1N.m (canonical, live copy)
% and chan_koop_yu2024_jbes_oisv/legacy/utility/getARh_approx1N.m (dead there).
% Function renamed getARh_approx1N -> init_approx1N.
% This function approximates the SV model using a linear Gaussian state 
% space model where the log chi^2 errors are modeled as N(-1.27,4.94)
%
% See:
% Chan, J.C.C. (2023). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, 235(2): 1419-1446.

function h_hat = init_approx1N(s2,muh,rhoh,sigh2)
    ystar = log(s2);
    T = length(s2);
    Hrhoh = speye(T) - rhoh*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
    HiSH = Hrhoh'*sparse(1:T,1:T,[1-rhoh^2, ones(1,T-1)])*Hrhoh;
    Kh = HiSH/sigh2 + speye(T)/4.94;
    h_hat = Kh\(muh/sigh2*HiSH*ones(T,1) + (ystar+1.27)/4.94);
end