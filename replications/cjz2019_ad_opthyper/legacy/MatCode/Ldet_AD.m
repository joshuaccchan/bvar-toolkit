% Support script for Chan, Jacobi and Zhu (2019)
%
% See: 
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019
function ld = Ldet_AD(Omega,k)
%ld=Sumvector(Log(DDiag(Cholasky2(Omega,K_nn,E_n,I_n,I_nn),k))); 
%ld.v=ld.v*2;
%ld.d=ld.d*2;
for i=1:size(Omega.d,2)
ld.d(1,i)=trace(Omega.inv*reshape(Omega.d(:,i),k,k));
end
ld.v=ldet(Omega.v);
end