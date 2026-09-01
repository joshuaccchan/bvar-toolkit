%This function compute the derivative of the determinant of a square matrix

function D=d_det(d,A,dA,n,k)
invA=A\speye(n);
 for i=1:k
 D(1,i)=d*trace(invA*reshape(dA(:,i),n,n));
end
