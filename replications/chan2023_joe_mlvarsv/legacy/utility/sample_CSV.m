% This script samples the common stochastic volatility h
%
% See:
% See:
% Chan, J.C.C. (2022). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, forthcoming.

function [h,is_accept] = sample_CSV(s2,rho,sigh2,h,n,is_ForcedAccept)
if nargin < 6
    is_ForcedAccept = false;
end
is_accept = 0;
T = size(s2,1);
Hrho = speye(T) - rho*sparse(2:T,1:(T-1),ones(1,T-1),T,T);
HiSH = Hrho'*sparse(1:T,1:T,[(1-rho^2)/sigh2; 1/sigh2*ones(T-1,1)])*Hrho;
errh = 1; ht = h;
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