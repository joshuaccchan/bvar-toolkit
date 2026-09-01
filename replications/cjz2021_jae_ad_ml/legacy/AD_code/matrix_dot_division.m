% This is written by Dan Zhu(dan.zhu@monash.edu) on the 16th of Jan 2020
% This function computes matrix element-wise dividsion as well as its derivatives.

function C=matrix_dot_division(A,B)
if isfield(A,'d')&&isfield(B,'d')==0
    [n,k]=size(B);d=n*k;
    C.v=A.v./B;
    C.d=sparse(1:d,1:d,1./B(:))*A.d;
    
elseif isfield(A,'d')==0&&isfield(B,'d')
    [n,k]=size(A);d=n*k;
    C.v=A./B.v;
    C.d=sparse(1:d,1:d,-C.v(:)./B.v(:))*B.d;
else
    [n,k]=size(A.v);d=n*k;
C.v=A.v./B.v;
C.d=sparse(1:d,1:d,1./B.v(:))*A.d+sparse(1:d,1:d,-C.v(:)./B.v)*B.d;
end