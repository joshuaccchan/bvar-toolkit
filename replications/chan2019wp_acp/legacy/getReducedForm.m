% This script converts the structrual-form draws to the reduced-form
% parameters
%
% See:
% Chan, J.C.C. (2019). Asymmetric conjugate priors for large Bayesian VARs,
% CAMA Working Papers 51/2019

function [store_Btilde,store_Sigtilde] = getReducedForm(store_alp,store_beta,store_Sig)
[nsim,n] = size(store_Sig);
k_beta = size(store_beta,2);
p = (k_beta/n-1)/n;
A_id = nonzeros(tril(reshape(1:n^2,n,n),-1)');
A = eye(n);
store_Btilde = zeros(nsim,n^2*p+n,1);
store_Sigtilde = zeros(nsim,n,n);
    % compute reduced-form parameters
for isim = 1:nsim
    alp = store_alp(isim,:)';
    beta = store_beta(isim,:)';
    sig = store_Sig(isim,:)';    
    
        % trasnform the parameters into reduced-form
    A(A_id) = alp;
    Sig = (A\sparse(1:n,1:n,sig))/A';    
    Btilde = (A\(reshape(beta,n*p+1,n)'))';

    store_Btilde(isim,:) = reshape(Btilde,n^2*p+n,1);
    store_Sigtilde(isim,:,:) = Sig;    
end
end

