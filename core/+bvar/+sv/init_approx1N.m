% bvar.sv.init_approx1N - crude 1-component log-chi2 approximation used to
% initialize SV paths.
%
%   h_hat = bvar.sv.init_approx1N(s2, muh, rhoh, sigh2)
%
%   s2                : T x 1 squared errors
%   muh, rhoh, sigh2  : mean, AR(1) coefficient and innovation variance of the
%                       log-volatility path
%   h_hat             : T x 1 conditional mean of h under the approximation
%
% The SV model is approximated by a linear Gaussian state space model in which
% the log chi^2 errors are modeled as N(-1.27,4.94), one normal in place of the
% 7-component mixture. Deterministic and deliberately crude: the result is a
% starting value for a sampler.
%
% Provenance and the legacy copies this stands in for: tests/variant_map.md.
% Equivalence: tests/unit/test_init_approx1N.m.
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