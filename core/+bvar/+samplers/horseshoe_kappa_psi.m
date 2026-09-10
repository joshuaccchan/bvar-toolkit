% bvar.samplers.horseshoe_kappa_psi - one sweep of the Minnesota-type HORSESHOE
% hierarchical shrinkage block: given the current coefficient vector theta
% (stacked k*n, intercept first per equation), draw the local scales psi
% (inverse-gamma), their auxiliaries z_psi, the global own-lag/other-lag scales
% kappa(1:2) (inverse-gamma), and their auxiliaries z_kappa, in exactly that
% order. kappa(3:4) pass through untouched (the OI model carries
% kappa(3) = NaN, the CS model kappa(3) = 1; kappa(4) = 100 intercepts in both).
%
%   [psi_kappa1,psi_kappa2,z_psi1,z_psi2,kappa,z_kappa] = ...
%       bvar.samplers.horseshoe_kappa_psi(theta,idx_kappa1,idx_kappa2,C, ...
%                                        kappa,z_psi1,z_psi2,z_kappa)
%
% Caller contract: the Psi reassembly, Psi(idx_kappa1) = psi_kappa1 and
% Psi(idx_kappa2) = psi_kappa2, stays with the caller. It consumes no rng and
% Psi is not read inside the block, so where the caller puts it in the sweep is
% draw-neutral.
%
% rng consumption, in order: gamrnd n*p-vector (psi1), gamrnd (n-1)*n*p-vector
% (psi2), gamrnd n*p-vector (z_psi1), gamrnd (n-1)*n*p-vector (z_psi2), two
% scalar gamrnd (kappa1, kappa2), one 1x2 gamrnd (z_kappa). NOTE the shape
% quirk: z_kappa may enter the FIRST call as the 2x1 column drawn at chain
% init, and leaves every call as a 1x2 row (gamrnd inherits the shape of
% kappa(1:2)).
%
% NEVER merge with bvar.samplers.gig_shrinkage: that is the MAHP normal-gamma
% (GIG) block, a different prior family with a different draw sequence.
% Provenance and the legacy copies this stands in for: tests/variant_map.md.
% Equivalence: tests/unit/test_oisv_equivalence.m.
%
% See:
% Chan, J.C.C., Koop, G. and Yu, X. (2024). Large Order-Invariant Bayesian
% VARs with Stochastic Volatility, Journal of Business and Economic
% Statistics, 42(2): 825-837.

function [psi_kappa1,psi_kappa2,z_psi1,z_psi2,kappa,z_kappa] = ...
    horseshoe_kappa_psi(theta,idx_kappa1,idx_kappa2,C,kappa,z_psi1,z_psi2,z_kappa)
np  = numel(idx_kappa1);    % = n*p
nnp = numel(idx_kappa2);    % = (n-1)*n*p

    % sample psi
tmpv1 = 1./z_psi1 + theta(idx_kappa1).^2./(2*C(idx_kappa1)*kappa(1));
tmpv2 = 1./z_psi2 + theta(idx_kappa2).^2./(2*C(idx_kappa2)*kappa(2));
psi_kappa1 = 1./gamrnd(1,1./tmpv1);
psi_kappa2 = 1./gamrnd(1,1./tmpv2);

    % sample z_psi
z_psi1 = 1./gamrnd(1,1./(1+1./psi_kappa1));
z_psi2 = 1./gamrnd(1,1./(1+1./psi_kappa2));

    % sample kappa1 and kappa2
tmpc1 = 1/z_kappa(1) + sum(theta(idx_kappa1).^2./(2*psi_kappa1.*C(idx_kappa1)));
tmpc2 = 1/z_kappa(2) + sum(theta(idx_kappa2).^2./(2*psi_kappa2.*C(idx_kappa2)));
kappa(1) = 1./gamrnd((np+1)/2,1./tmpc1);
kappa(2) = 1./gamrnd((nnp+1)/2,1./tmpc2);

    % sample z_kappa
z_kappa = 1./gamrnd(1,1./(1+1./kappa(1:2)));
end
