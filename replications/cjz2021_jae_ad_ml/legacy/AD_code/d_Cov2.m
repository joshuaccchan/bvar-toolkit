% This is written by Dan Zhu(dan.zhu@monash.edu) on the 6th of Jan 2020
% This function computes covariance matrix of a sample as well as its derivatives.

function C=d_Cov2(X)
[n,k]=size(X.v);
m.v=mean(X.v);
m.d=mean(X.d,3);
error.v=X.v;
error.d=zeros(k*n,size(m.d,2));
for i=1:n
  ind=i+(0:k-1)*n;
   error.d(ind,:)=X.d(:,:,i);
end    
error=matrix_minus_vector(error,m,1);
C=Mtimes(Mtrans(error,[]),error);
C=MtimeScalar(C,1/n);



