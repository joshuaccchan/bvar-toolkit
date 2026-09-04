% bvar.sv.sv0_params - posterior draw of the SV state-equation parameters
% (phi, sig2) for stationary zero-mean AR(1) log-volatilities.
%
% Extracted 2026-09-01 (step 4, SV/prior core). Canonical source (body verbatim):
% chan_koop_yu2024_jbes_oisv/legacy/utility/sample_SV0para.m (single copy).
% Edits made: function renamed sample_SV0para -> sv0_params; the phi-candidate
% MH truncation bound (hard-coded .99 in the canonical copy) is promoted to the
% optional 4th argument phi_bnd, DEFAULTING to .99, so the default call
% reproduces the OISV copy exactly - same outputs and the same sequence and
% count of gamrnd/randn/rand calls under the same rng seed.
%
% This is NOT a special case of bvar.sv.sv_params to be merged away: the OISV
% pair deliberately splits the zero-mean sampler into its own file with a
% DIFFERENT phi truncation bound (.99 here vs .999 in sample_SVpara); see the
% never-merge section of tests/variant_map.md.
%
% This function samples the SV parameters phi, and sig2

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
