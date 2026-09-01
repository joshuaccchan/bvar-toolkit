% This function samples the parameters in the state equation for the 
% log-volatilities
%
% See:
% Chan, J.C.C. (2022). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, forthcoming.

function [mu,phi,sig2,flag_phi] = sample_SVpara(h,mu,phi,Hyper)    
[T,n] = size(h);
    % sample sig2
tmp_h = h-repmat(mu',T,1);
e_h = [tmp_h(1,:).*sqrt(1-phi.^2)'; tmp_h(2:end,:)-repmat(phi',T-1,1).*tmp_h(1:end-1,:)];
sig2 = 1./gamrnd(Hyper.nuh+T/2,1./(Hyper.Sh + sum(e_h.^2)'/2));

    % sample phi
Kphi = 1./Hyper.Vphi + sum(tmp_h(1:T-1,:).^2)'./sig2;
phi_hat = (Hyper.phi0./Hyper.Vphi + sum(tmp_h(1:T-1,:).*tmp_h(2:T,:))'./sig2)./Kphi;
phic = phi_hat + 1./sqrt(Kphi).*randn(n,1);
flag_phi = zeros(n,1);
for ii = 1:n
    g_phi = @(x) .5*log(1-x^2) -.5*(1-x^2)/sig2(ii)*tmp_h(1,ii)^2;
    if abs(phic(ii))<.998
        alpMH = exp(g_phi(phic(ii))-g_phi(phi(ii)));
        if alpMH>rand
            phi(ii) = phic(ii);
            flag_phi(ii) = 1;
        end
    end 
end   

if mu ~= 0  % sample mu only if it is not assumed to be zero
        % sample mu    
    Kmu = 1./Hyper.Vmu + ((1-phi.^2) + (T-1)*(1-phi).^2)./sig2;
    mu_hat = (Hyper.mu0./Hyper.Vmu + (1-phi.^2)./sig2.*h(1,:)' + ...
       (1-phi)./sig2.*sum(h(2:end,:)-repmat(phi',T-1,1).*h(1:end-1,:))')./Kmu;
    mu = mu_hat + 1./sqrt(Kmu).*randn(n,1);
end
end