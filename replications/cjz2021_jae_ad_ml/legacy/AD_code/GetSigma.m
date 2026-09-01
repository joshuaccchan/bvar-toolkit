function [mu,Sigma_0]=GetSigma(m, v, k)
global L
mu.v=zeros(m,1);
mu.d=zeros(m,L);
Sigma_0.v=1/v*eye(m);
Sigma_0.d=zeros(m,L);
Sigma_0.d(:,k)=1/v;
Sigma_0.d=d_diag(Sigma_0.d,m); 
Sigma_0=dmatrix(Sigma_0,m,L);
end
       

