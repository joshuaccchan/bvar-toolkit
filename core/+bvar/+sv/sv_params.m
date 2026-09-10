% bvar.sv.sv_params - posterior draw of the SV state-equation parameters
% (mu, phi, sig2) for stationary AR(1) log-volatilities with mean mu.
%
%   [mu,phi,sig2,flag_phi] = bvar.sv.sv_params(h, mu, phi, Hyper)
%   [mu,phi,sig2,flag_phi] = bvar.sv.sv_params(h, mu, phi, Hyper, phi_bnd)
%
%   h        : T x (n+r) log-volatility paths. The first n columns have mean mu
%              (n = numel(mu)); the last r columns are zero-mean and share the
%              same phi and sig2 draws
%   mu       : n x 1 current means; the new draw on output
%   phi      : (n+r) x 1 current AR(1) coefficients; the new draw on output
%   Hyper    : struct with nuh, Sh (inverse-gamma prior on sig2), phi0, Vphi
%              (normal prior on phi) and mu0, Vmu (normal prior on mu)
%   phi_bnd  : truncation bound on the phi candidate, which can be accepted
%              only if |phic| < phi_bnd; default .999, the OISV canonical
%              value, while ml_varsv uses .998
%   sig2     : (n+r) x 1 innovation variances
%   flag_phi : (n+r) x 1, 1 where the phi candidate was accepted
%
% Trap: the mu block is gated on `if mu~=0`, which under MATLAB semantics runs
% only when EVERY element of mu is nonzero; a single exact zero silently skips
% the mu draw for all n series. Callers initialize mu from the data, so it is
% all-nonzero in practice and mu is drawn every sweep.
%
% NEVER merge with bvar.sv.sv0_params, the separate zero-mean sampler, which
% carries a different truncation bound (.99).
%
% rng consumption: one gamrnd for sig2, one randn(n+r,1) for the phi candidates,
% then one rand per candidate falling inside phi_bnd - so that count is
% data-dependent - and finally randn(n,1) for mu when the gate above passes.
%
% Provenance and the legacy copies this stands in for: tests/variant_map.md.
% Equivalence: tests/unit/test_sv_params.m, test_sv_params_mlvarsv.m.
%
% See:
% Chan, J.C.C., Koop, G. and Yu, X. (2024). Large Order-Invariant Bayesian
% VARs with Stochastic Volatility, Journal of Business and Economic
% Statistics, 42(2): 825-837.

function [mu,phi,sig2,flag_phi] = sv_params(h,mu,phi,Hyper,phi_bnd)
if nargin < 5
    phi_bnd = .999;     % OISV canonical truncation bound
end
n = size(mu,1);
[T,npr] = size(h);
r = npr - n;

    % sample sig2
tmp_h = [h(:,1:n)-repmat(mu',T,1), h(:,n+1:end)];
e_h = [tmp_h(1,:).*sqrt(1-phi.^2)'; tmp_h(2:end,:)-repmat(phi',T-1,1).*tmp_h(1:end-1,:)];
sig2 = 1./gamrnd(Hyper.nuh+T/2,1./(Hyper.Sh + sum(e_h.^2)'/2));

    % sample phi
Kphi = 1./Hyper.Vphi + sum(tmp_h(1:T-1,:).^2)'./sig2;
phi_hat = (Hyper.phi0./Hyper.Vphi + sum(tmp_h(1:T-1,:).*tmp_h(2:T,:))'./sig2)./Kphi;
phic = phi_hat + 1./sqrt(Kphi).*randn(n+r,1);
flag_phi = zeros(n+r,1);
for ii = 1:n+r
    g_phi = @(x) .5*log(1-x^2) -.5*(1-x^2)/sig2(ii)*tmp_h(1,ii)^2;
    if abs(phic(ii))<phi_bnd
        alpMH = exp(g_phi(phic(ii))-g_phi(phi(ii)));
        if alpMH>rand
            phi(ii) = phic(ii);
            flag_phi(ii) = 1;
        end
    end
end

    % sample mu
if mu~=0
Kmu = 1./Hyper.Vmu + ((1-phi(1:n).^2) + (T-1)*(1-phi(1:n)).^2)./sig2(1:n);
mu_hat = (Hyper.mu0./Hyper.Vmu + (1-phi(1:n).^2)./sig2(1:n).*h(1,1:n)' ...
   +(1-phi(1:n))./sig2(1:n).*sum(h(2:end,1:n)-repmat(phi(1:n)',T-1,1).*h(1:end-1,1:n))')./Kmu;
mu = mu_hat + 1./sqrt(Kmu).*randn(n,1);
end
end
