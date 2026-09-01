% This function computes the derivative of the determinant of a square matrix
%
% See: 
% Chan, J. C. C., L. Jacobi, and D. Zhu (2019). Efficient Selection of
% Hyperparameters in Large Bayesian VARs Using Automatic Differentiation, 
% CAMA Working Paper 46/2019

function D=d_det(d,A,dA,n,k)
invA=A\speye(n);
for i=1:k
D(1,i)=d*trace(invA*reshape(dA(:,i),n,n));
end