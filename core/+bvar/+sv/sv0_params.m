% bvar.sv.sv0_params - posterior draw of the SV state-equation parameters
% (phi, sig2) for stationary zero-mean AR(1) log-volatilities.
%
%   [phi,sig2,flag_phi] = bvar.sv.sv0_params(h, phi, Hyper)
%   [phi,sig2,flag_phi] = bvar.sv.sv0_params(h, phi, Hyper, phi_bnd)
%
%   h        : T x n zero-mean log-volatility paths, one column per series
%   phi      : n x 1 current AR(1) coefficients; the new draw on output
%   Hyper    : struct with nuh, Sh (inverse-gamma prior on sig2) and phi0,
%              Vphi (normal prior on phi)
%   phi_bnd  : truncation bound on the phi candidate, which can be accepted
%              only if |phic| < phi_bnd; default .99, the OISV canonical value
%   sig2     : n x 1 innovation variances
%   flag_phi : n x 1, 1 where the phi candidate was accepted
%
% NEVER merge with bvar.sv.sv_params. This is not that sampler at mu = 0: the
% OISV pair keeps the zero-mean sampler separate, with a different truncation
% bound (.99 here, .999 there).
%
% See:
% Chan, J.C.C., Koop, G. and Yu, X. (2024). Large Order-Invariant Bayesian
% VARs with Stochastic Volatility, Journal of Business and Economic
% Statistics, 42(2): 825-837.

function [phi,sig2,flag_phi] = sv0_params(h,phi,Hyper,phi_bnd)
if nargin < 4
    phi_bnd = .99;      % OISV canonical truncation bound
end
[T,n] = size(h);

    % sample sig2
e_h = [h(1,:).*sqrt(1-phi.^2)'; h(2:end,:)-repmat(phi',T-1,1).*h(1:end-1,:)];
sig2 = 1./gamrnd(Hyper.nuh+T/2,1./(Hyper.Sh + sum(e_h.^2)'/2));

    % sample phi
Kphi = 1./Hyper.Vphi + sum(h(1:T-1,:).^2)'./sig2;
phi_hat = (Hyper.phi0./Hyper.Vphi + sum(h(1:T-1,:).*h(2:T,:))'./sig2)./Kphi;
phic = phi_hat + 1./sqrt(Kphi).*randn(n,1);
flag_phi = zeros(n,1);
for ii = 1:n
    g_phi = @(x) .5*log(1-x^2) -.5*(1-x^2)/sig2(ii)*h(1,ii)^2;
    if abs(phic(ii))<phi_bnd
        alpMH = exp(g_phi(phic(ii))-g_phi(phi(ii)));
        if alpMH>rand
            phi(ii) = phic(ii);
            flag_phi(ii) = 1;
        end
    end
end
end
