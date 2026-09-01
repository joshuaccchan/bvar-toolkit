% This is written by Dan Zhu(dan.zhu@monash.edu) on the 6th of June 2019
% This function computes matrix minus a vector, 
% if b is rwo vector, C(i,:)=A(i,:)-b;
% if b is column vector  C(:,i)=A(:,i)-b;
%as well as its derivatives.

function C=matrix_minus_vector(A,b,row)
if row==1
if isfield(A,'d')&& isfield(b,'d')==0
    n=size(A.v,1);
    C.v=A.v-ones(n,1)*b;
    C.d=A.d;
elseif isfield(A,'d')==0&& isfield(b,'d')
    [n,k]=size(A);
    C.v=A-ones(n,1)*b.v;
    C.d=-kron(speye(k),ones(n,1))*b.d;
else
   [n,k]=size(A.v);
C.v=A.v-ones(n,1)*b.v;
C.d=A.d-kron(speye(k),ones(n,1))*b.d;
end

else
    
if isfield(A,'d')&& isfield(b,'d')==0
    d=size(A.v,2);
    C.v=A.v-b*ones(1,d);
    C.d=A.d;
elseif isfield(A,'d')==0&& isfield(b,'d')
    [n,d]=size(A);
    C.v=A-b.v*ones(1,d);
    C.d=-kron(ones(d,1),speye(n))*b.d;
else
    [n,d]=size(A.v);
C.v=A.v-b.d*ones(1,d);
C.d=A.d-kron(ones(d,1),speye(n))*b.d;
end   
    
end    
 



