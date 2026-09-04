% bvar.sv.csv_armh - one AR-MH (accept-reject Metropolis-Hastings) sweep for the
% common stochastic volatility (CSV) log-volatility path h, marginal of the mixture
% indicators: Newton-Raphson mode-finding, Gaussian AR proposal, then MH correction.
%
% Extracted 2026-09-01 (step 4, SV/prior core). Canonical source:
%   chan2023_joe_mlvarsv/legacy/utility/sample_CSV.m  (body verbatim; renamed
%   sample_CSV -> csv_armh).
% Also canonicalizes exactly (verified draw-for-draw under fixed seeds, identical
% rand/randn call sequence and count, identical terminal rng state):
%   chan2020_jbes_kronecker/legacy/sample_h.m
%   chan2020_jbes_kronecker/legacy/realtime_forecasts/sample_h.m
%   chan2020_springer_largebvar/legacy/sample_h.m
% The three sample_h copies have executably identical bodies (differences are
% comments, `[h is_accept]` output-list syntax, a stray semicolon after `while`,
% and `10^(-3)` vs `1e-3` - bitwise-equal doubles). sample_CSV differs only by the
% is_ForcedAccept flag; with the flag false (the default here) `log(rand)` is still
% evaluated before the short-circuit `||`, so the rng sequence matches sample_h
% exactly. Parameterized: is_ForcedAccept (optional, default false = sample_h
% behavior; true = always move to the proposal, as ml_varsv's VAR_CSV.m does at
% initialization and during early burn-in iterations).
%
% Second parameterization added 2026-09-02 (step 8, Kronecker family pass):
% ht_start (optional, default h = the legacy sample_h/sample_CSV behavior) -
% the Newton-Raphson mode-search starting point. The inline h step of the
% marginal-likelihood reduced run in chan2020_jbes_kronecker/legacy/
% ml_BVAR_CSV.m (lines 51-88) is byte-equivalent to this body EXCEPT that it
% starts the NR search at the posterior mean path h_mean instead of the
% current h (its line 52: ht = h_mean) and force-accepts on the reduced run's
% first sweep; csv_armh(s2,rho,sigh2,h,n,isim==1,h_mean) reproduces it
% bitwise, including the rng call sequence (log(rand) is evaluated before the
% short-circuit ||, exactly as in the legacy `if alpMH > log(rand) || isim == 1`).
% Default calls are unchanged bit-for-bit.
%
% Inputs:  s2    - T x 1, sum over the n series of squared (orthogonalized,
%                  lambda-scaled) errors at each t
%          rho   - AR(1) coefficient of h
%          sigh2 - innovation variance of h
%          h     - T x 1 current log-volatility path
%          n     - number of series
% Outputs: h, is_accept (1 if the MH proposal was accepted)
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error
% covariance structure, Journal of Business and Economic Statistics, 38(1), 68-79.
% Chan, J.C.C. (2023). Comparing stochastic volatility specifications for large
% Bayesian VARs, Journal of Econometrics, 235(2), 1419-1446.

function [h,is_accept] = csv_armh(s2,rho,sigh2,h,n,is_ForcedAccept,ht_start)
if nargin < 6
    is_ForcedAccept = false;
end
if nargin < 7
    ht_start = h;
end
is_accept = 0;
T = size(s2,1);
Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
HiSH = Hrho'*sparse(1:T,1:T,[(1-rho^2)/sigh2; 1/sigh2*ones(T-1,1)])*Hrho;
errh = 1; ht = ht_start;
while errh> 1e-3
    eht = exp(ht);
    sieht = s2./eht;
    fh = -n/2 + .5*sieht;
    Gh = .5*sieht;
    Kh = HiSH + sparse(1:T,1:T,Gh);
    newht = Kh\(fh+Gh.*ht);
    errh = max(abs(newht-ht));
    ht = newht;
end
CKh = chol(Kh,'lower');
% AR-step
hstar = ht;
logc = -.5*hstar'*HiSH*hstar - n/2*sum(hstar) - .5*exp(-hstar)'*s2 + log(3);
flag = 0;
while flag == 0
    hc = ht + CKh'\randn(T,1);
    alpARc = -.5*hc'*HiSH*hc - n/2*sum(hc) - .5*exp(-hc)'*s2 ...
        + .5*(hc-ht)'*Kh*(hc-ht) - logc;
    if alpARc > log(rand)
        flag = 1;
    end
end
% MH-step
alpAR = -.5*h'*HiSH*h - n/2*sum(h) -.5*exp(-h)'*s2 + .5*(h-ht)'*Kh*(h-ht) - logc;
if alpAR < 0
    alpMH = 1;
elseif alpARc < 0
    alpMH = - alpAR;
else
    alpMH = alpARc - alpAR;
end
if alpMH > log(rand) || is_ForcedAccept
    h = hc;
    is_accept = 1;
end
end
