function draw = anormrnd(mu,rho)
w = 1/(1+exp(2*mu/rho));
if w > rand
    mu1 = mu/2-sqrt(mu^2+4)/2;    
    sig21 = mu1^2*rho/(1+mu1^2);
    draw = mu1 + sqrt(sig21)*randn;
else
    mu2 = mu/2+sqrt(mu^2+4)/2;
    sig22 = mu2^2*rho/(1+mu2^2);
    draw = mu2 + sqrt(sig22)*randn;
end

end