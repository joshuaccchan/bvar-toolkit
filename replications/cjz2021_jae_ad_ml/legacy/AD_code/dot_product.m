% This is written by Dan Zhu(dan.zhu@monash.edu) on the 16th of June 2018
% It computes derivative of the operation dot_prod of two column vector

function C=dot_product(A,B)
if isfield(A,'d')&& isfield(B,'d')==0
    C.v=dot(A.v,B);
    C.d=B'*A.d;
elseif isfield(A,'d')==0&& isfield(B,'d')
    C.v=dot(A,B.v);
    C.d=A'*B.d;
else
C.v=dot(A.v,B.v);
C.d=B.v'*A.d+A.v'*B.d;
end