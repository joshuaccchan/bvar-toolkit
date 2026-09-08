% This function obtains the posterior draws of B and Sigma
% under the natural conjugate prior

function [store_B,store_Sig] = sample_BSig_NCP(Y0,Y,p,prior,nsim)
[T,n] = size(Y);
tmpY = [Y0(end-p+1:end,:); Y];
Z = zeros(T,n*p); 
for ii=1:p
    Z(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
end
Z = [ones(T,1) Z];
k = n*p+1;
    % compute the parameters for the posterior density
ZZ = Z'*Z;
Btilde = ZZ\(Z'*Y);
KB = sparse(1:k,1:k,1./prior.VB) + ZZ;
    % posterior mean of the VAR coefficients, arranged as a k by n matrix
Bhat = KB\(sparse(1:k,1:k,prior.VB)\prior.B0 + ZZ*Btilde); 
Shat = prior.S0 + prior.B0'*sparse(1:k,1:k,1./prior.VB)*prior.B0 + Y'*Y - Bhat'*KB*Bhat;
Shat = (Shat+Shat')/2;
CKB = chol(KB,'lower');
store_B = zeros(nsim,n*k);
store_Sig = zeros(nsim,n,n);
for isim = 1:nsim
    Sig = iwishrnd(Shat,prior.nu0+T); 
    CSig = chol(Sig,'lower');
    B = Bhat + (CKB'\randn(k,n))*CSig';    
    store_B(isim,:) = B(:);
    store_Sig(isim,:,:) = Sig;
end   
end