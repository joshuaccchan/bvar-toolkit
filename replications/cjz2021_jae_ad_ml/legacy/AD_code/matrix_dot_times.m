% This is written by Dan Zhu(dan.zhu@monash.edu) on the 16th of Jan 2020
% This function computes matrix element-wise product as well as its derivatives.
function C=matrix_dot_times(A,B)
if isfield(A,'d')&&isfield(B,'d')==0
    [n,k]=size(B);d=n*k;
    C.v=A.v.*B;
    C.d=sparse(1:d,1:d,B(:))*A.d;
    
elseif isfield(A,'d')==0&&isfield(B,'d')
    [n,k]=size(A);d=n*k;
    C.v=A.*B.v;
    C.d=sparse(1:d,1:d,A(:))*B.d;
else
    [n,k]=size(A.v);d=n*k;
C.v=A.v.*B.v;
C.d=sparse(1:d,1:d,B.v(:))*A.d+sparse(1:d,1:d,A.v(:))*B.d;
end



