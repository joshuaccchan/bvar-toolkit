%This is written by Dan Zhu(dan.zhu@monash.edu) on the 4th of Feb that
%computes the vector sum of matrix over the row or column and its associate
%derivatives 

function B=Sum_matrix(A,dim)
[n,k]=size(A.v);

if dim==2
   B.v=sum(A.v,2);
   B.d=repmat(speye(n),1,k)*A.d;   
else
   B.v=sum(A.v);
   B.d=kron(speye(k),ones(1,n))*A.d;
end