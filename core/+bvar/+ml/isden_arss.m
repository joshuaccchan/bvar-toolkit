% bvar.ml.isden_arss - Gaussian importance-sampling density for one log-volatility
% path, fitted to its posterior draws under the state space parameterization
%   h_t = a_t + rho h_{t-1} + u_t,  u_t ~ N(0,b_t).
% rho is the maximizer of the concentrated likelihood over (-.99,.99); the
% per-period a_t, b_t are the mean and variance of the rho-differenced draws.
% Returns the mean h_hat and the PRECISION Kh_hat (band, sparse).
%
% Extracted 2026-09-03 (step 10, ml_varsv marginal-likelihood pass). Canonical
% source: chan2023_joe_mlvarsv/legacy/utility/getISden_ARSS.m (body and its
% concen_like_h subfunction verbatim; the only copy in the repo - no other
% package ships it). Called by all four bvar.ml.mlvarsv_* routines, once per
% series; no rng. Verified bitwise by tests/unit/test_mlvarsv_ml_densities.m
% and end-to-end by tests/unit/test_mlvarsv_ml.m.
%
%   [h_hat,Kh_hat,r_hat,m_hat,v_hat] = bvar.ml.isden_arss(store_h)
%   store_h - R x T matrix of posterior draws (one row per draw)
%
% See:
% Chan, J.C.C. (2023). Comparing stochastic volatility specifications for large
% Bayesian VARs, Journal of Econometrics, 235(2), 1419-1446.

function [h_hat,Kh_hat,r_hat,m_hat,v_hat] = isden_arss(store_h)
T = size(store_h,2);
r_hat = fminbnd(@(x) -concen_like_h(store_h,x),-.99,.99);
Hr = speye(T) -r_hat*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
[~,m_hat,v_hat] = concen_like_h(store_h,r_hat);
h_hat = Hr\m_hat;
Kh_hat = Hr'*sparse(1:T,1:T,1./v_hat)*Hr;

end

function [concen_like,m_hat,v_hat] = concen_like_h(store_h,rho)
[R,T] = size(store_h);
Hrho = speye(T) -rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
Hrhoh = Hrho*store_h';
m_hat = mean(Hrhoh,2);
E = Hrhoh - repmat(m_hat,1,R);
v_hat = mean(E.^2,2);

concen_like = -T*R/2*log(2*pi) -R/2*sum(log(v_hat)) ...
    -.5*sum(sum(E.^2./repmat(v_hat,1,R)));
end
