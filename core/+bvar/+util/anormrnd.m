% bvar.util.anormrnd - one draw from the two-component (bimodal) normal
% approximation used for the FIRST rotation coordinate zeta(1) in the row-wise
% B0 step of the order-invariant SVAR-SV sampler: mixture weight
% w = 1/(1+exp(2*mu/rho)) on the negative mode mu1 = mu/2 - sqrt(mu^2+4)/2
% (else the positive mode mu2 = mu/2 + sqrt(mu^2+4)/2), each with variance
% muj^2*rho/(1+muj^2).
%
%   draw = bvar.util.anormrnd(mu, rho)
%
%   mu, rho : scalar parameters of the approximation above
%   draw    : scalar draw
%
% rng consumption: exactly one rand THEN one randn per call.
%
% Not a truncated normal - do NOT fold into bvar.util.tnormrnd (different
% density, different rng sequence).
%
% Provenance and the legacy copies this stands in for: tests/variant_map.md.
% Equivalence: tests/unit/test_anormrnd.m.
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
