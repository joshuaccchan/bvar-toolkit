% This is written by Dan Zhu(dan.zhu@monash.edu) on the 16th of Jan 2020
% This function computes matrix multiplication as well as its derivatives.
function [C] = Mtimes(A, B)
if isfield(A,'d')&&isfield(B,'d')==0
    C.d=sparse(kron(B', speye(size(A.v,1)))) * A.d;
    C.v=A.v*B;
elseif isfield(A,'d')==0&&isfield(B,'d')
    C.d=sparse( kron(speye(size(B.v,2)), A)) * B.d;
    C.v=A*B.v;
else  
  I_a = speye(size(A.v,1));
  I_b =speye(size(B.v,2));
  f1 =sparse( kron(I_b, A.v));
  f2 = sparse(kron(B.v', I_a));
  C.d = f1 * B.d + f2 * A.d;
  C.v=A.v*B.v;
end
