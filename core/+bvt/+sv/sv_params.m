% bvt.sv.sv_params - posterior draw of the SV state-equation parameters
% (mu, phi, sig2) for stationary AR(1) log-volatilities with mean mu.
%
% Extracted 2026-09-01 (step 4, SV/prior core). Canonical source (body verbatim):
% chan_koop_yu2024_jbes_oisv/legacy/utility/sample_SVpara.m (the only OISV copy).
% Edits made: function renamed sample_SVpara -> sv_params; the phi-candidate MH
% truncation bound (hard-coded .999 in the canonical copy) is promoted to the
% optional 5th argument phi_bnd, DEFAULTING to .999, so the default call
% reproduces the OISV copy exactly - same outputs and the same sequence and
% count of gamrnd/randn/rand calls under the same rng seed.
%
% NOT canonicalized here (never-merge, see tests/variant_map.md):
% chan2023_joe_mlvarsv/legacy/utility/sample_SVpara.m - same name but
% numerically different: phi bound .998, no n+r column split (mu is subtracted
% from ALL columns of h), and full-vector mu/phi/sig2 in the mu block.
%
% Notes on the verbatim body:
% - h is T x (n+r): the first n columns have mean mu (n = numel(mu)); the last
%   r columns are zero-mean log-volatilities sharing the same phi/sig2 draws.
% - `if mu~=0` on a vector follows MATLAB semantics: the mu block runs only if
%   EVERY element of mu is nonzero; a single exact zero skips the whole mu draw.
%   (Header corrected 2026-09-02, step 7: CS_MH initializes mu(ii) =
%   mean(log(s2i)) from the data - lines 34-37 - so mu is almost surely
%   all-nonzero, the gate passes, and mu IS updated every sweep; verified
%   draw-for-draw by tests/unit/test_oisv_equivalence.m, whose store_hpara
%   comparison includes the mu draws. The gate is kept verbatim regardless.)
%
% This function samples the SV parameters mu, phi, and sig2

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
