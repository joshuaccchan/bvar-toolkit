% bvar.util.anormrnd - one draw from the two-component (bimodal) normal
% approximation used for the FIRST rotation coordinate zeta(1) in the row-wise
% B0 step of the order-invariant SVAR-SV sampler: mixture weight
% w = 1/(1+exp(2*mu/rho)) on the negative mode mu1 = mu/2 - sqrt(mu^2+4)/2
% (else the positive mode mu2 = mu/2 + sqrt(mu^2+4)/2), each with variance
% muj^2*rho/(1+muj^2). Consumes exactly one rand THEN one randn per call.
%
% Extracted 2026-09-02 (step 7, OISV family pass). Canonical source (body
% verbatim): chan_koop_yu2024_jbes_oisv/legacy/utility/anormrnd.m (single copy;
% called by SVARSV_MH.m line 63 and forecast_SVARSV_MH.m line 57, both via
% bvar.structural.b0_row_sampler's canonical block). Only the namespace was
% added; nothing renamed or parameterized. Model-specific draw - do NOT fold
% into bvar.util.tnormrnd (a truncated-normal sampler; entirely different
% density and rng sequence).
%
% See:
% Chan, J.C.C., Koop, G. and Yu, X. (2024). Large Order-Invariant Bayesian
% VARs with Stochastic Volatility, Journal of Business and Economic
% Statistics, 42(2): 825-837.

function draw = anormrnd(mu,rho)
w = 1/(1+exp(2*mu/rho));
if w > rand
    mu1 = mu/2-sqrt(mu^2+4)/2;
    sig21 = mu1^2*rho/(1+mu1^2);
    draw = mu1 + sqrt(sig21)*randn;
else
    mu2 = mu/2+sqrt(mu^2+4)/2;
    sig22 = mu2^2*rho/(1+mu2^2);
    draw = mu2 + sqrt(sig22)*randn;
end

end
