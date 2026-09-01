% This script samples the log-volatility (random walk state equation) 

function h = sample_SVRW(ystar,h,sig2,h0)
T = length(h);
    % 7-component normal mixture
p_N = [0.0073 .10556 .00002 .04395 .34001 .24566 .2575];
m_N = [-10.12999 -3.97281 -8.56686 2.77786 .61942 1.79518 -1.08819] - 1.2704;  %% means already adjusted!! %%
sig2_N = [5.79596 2.61369 5.17950 .16735 .64009 .34023 1.26261];

    % sample S from a 7-point distrete distribution
tmprand = rand(T,1);
q = repmat(p_N,T,1).*normpdf(repmat(ystar,1,7),repmat(h,1,7)+repmat(m_N,T,1),...
    repmat(sqrt(sig2_N),T,1));
q = q./repmat(sum(q,2),1,7);
S = 7 - sum(repmat(tmprand,1,7)<cumsum(q,2),2)+1;
    
    % sample h
% y^* = h + d + \epsilon, \epsilon \sim N(0,\Omega),
% Hh = \alpha + \nu, \nu \ sim N(0,S),
% where d_t = Ez_t, \Omega = diag(\omega_1,...,\omega_n), 
% \omega_t = var z_t, S = diag(sig, \ldots, sig)

H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);
HiSH = H'*sparse(1:T,1:T,1/sig2*ones(T,1))*H;
d = m_N(S)'; iOmega = spdiags(1./sig2_N(S)',0,T,T);
Kh = HiSH + iOmega;
muh = H\[h0;sparse(T-1,1)];  
hhat = Kh\(HiSH*muh + iOmega*(ystar-d));
h = hhat + chol(Kh,'lower')'\randn(T,1);
end
