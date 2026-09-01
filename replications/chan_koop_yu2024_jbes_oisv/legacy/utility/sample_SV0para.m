% This function samples the SV parameters phi, and sig2

function [phi,sig2,flag_phi] = sample_SV0para(h,phi,Hyper)
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
    if abs(phic(ii))<.99
        alpMH = exp(g_phi(phic(ii))-g_phi(phi(ii)));
        if alpMH>rand
            phi(ii) = phic(ii);
            flag_phi(ii) = 1;
        end
    end 
end
end