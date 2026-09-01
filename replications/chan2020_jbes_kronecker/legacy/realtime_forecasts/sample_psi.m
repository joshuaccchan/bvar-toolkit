function [psi,flag,psi_hat,Kpsic] = sample_psi(psi,e,sig2,lpri_psi)

lp_psi = @(x) llike_uMA1(x,e,sig2) + lpri_psi(x);
psi_hat = fminbnd(@(x)-lp_psi(x),-.99,.99);

    % compute Hessian at psi_hat
Del = .01; % step size
f_1 = llike_uMA1(psi_hat+Del,e,sig2);
f_0 = llike_uMA1(psi_hat,e,sig2);
f_n1 = llike_uMA1(psi_hat-Del,e,sig2);
d2f_psi = (f_n1 - 2*f_0 + f_1)/Del^2;
Kpsic = - d2f_psi;        

    % sample psi
[Cpsi, flag] = chol(Kpsic,'lower');
if flag ~= 0
    Cpsi = .1;
end    
psic = psi_hat + Cpsi'\randn; 
if abs(psic) < .99
    lg = @(x) -.5*(x-psi_hat)'*Kpsic*(x-psi_hat); % kernel of the proposal
    alpMH = llike_uMA1(psic,e,sig2) - llike_uMA1(psi,e,sig2) ...
        + lpri_psi(psic) - lpri_psi(psi) + lg(psi) - lg(psic); 
else
    alpMH = -inf;
end
flag = alpMH > log(rand);
if flag
    psi = psic;
end 
end

