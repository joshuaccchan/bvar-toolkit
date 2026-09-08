% This function computes the imuplse responses
%
% See:
% Chan, J.C.C. (2022). Asymmetric Conjugate Priors for Large Bayesian VARs,
% Quantitative Economics, 13(3): 1145-1169
%
% Inputs:
% A: reduced-form coef., each column contains coef. for each equation
% L: L*L' = reduced-form covariance matrix
%
% Output:
% response(:,:,it) contains the impulse responses at horizon it-1

function response = IRredu(A,L,nstep,nshock)
[np, n] = size(A);
p = np/n;
response = zeros(n,nshock,nstep);
Acomp = [A'; sparse(1:n*(p-1),1:n*(p-1),ones(n*(p-1),1),n*(p-1),np)];
response(:,:,1) = L(:,1:nshock);
Apower = Acomp;
for it = 2:nstep
    response(:,:,it) = Apower(1:n,1:n)*L(:,1:nshock);
    Apower = Apower*Acomp;
end
end