function B=Mstack(A,n,m, I1,I2)
B.v=sparse(n,m);
B.v(I1)=A.v(I2);
B.d=sparse(n*m,size(A.d,2));
B.d(I1,:)=A.d(I2,:);