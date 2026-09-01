% This script samples the log-volatilities (AR state equation) 

function [h,S] = sample_SV(ystar,h,mu,rho,sig2)
T = length(h);
    % 7-component normal mixture
p_N = [0.0073 .10556 .00002 .04395 .34001 .24566 .2575];
m_N = [-10.12999 -3.97281 -8.56686 2.77786 .61942 1.79518 -1.08819] - 1.2704;  % means already adjusted!!
sig2_N = [5.79596 2.61369 5.17950 .16735 .64009 .34023 1.26261];

    % sample S from a 7-point distrete distribution
tmprand = rand(T,1);
q = repmat(p_N,T,1).*normpdf(repmat(ystar,1,7),repmat(h,1,7)+repmat(m_N,T,1),...
    repmat(sqrt(sig2_N),T,1));
q = q./repmat(sum(q,2),1,7);
S = 7 - sum(repmat(tmprand,1,7)<cumsum(q,2),2)+1;
    
    % sample h using the precision sampler
% y^* = h + d + \epsilon, \epsilon \sim N(0,\Omega)
% h  ~ N(mu 1_T, sig2(Hrho' S^{-1} Hrho)^{-1}),
% where S^{-1} = diag(1-rho^2,1,1,...,1)

Hrho = speye(T) - sparse(2:T,1:(T-1),rho*ones(1,T-1),T,T);    
HiSH = Hrho'*sparse(1:T,1:T,[1-rho^2, ones(1,T-1)])*Hrho;
d = m_N(S)'; iOmega = sparse(1:T,1:T,1./sig2_N(S));
Kh = HiSH/sig2 + iOmega;
h_hat = Kh\(mu/sig2*HiSH*ones(T,1) + iOmega*(ystar-d));
h = h_hat + chol(Kh,'lower')'\randn(T,1);
end
