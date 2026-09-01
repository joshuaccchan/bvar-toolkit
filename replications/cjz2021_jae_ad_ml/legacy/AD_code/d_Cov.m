function D=d_Cov(X,dX)
m=mean(X);
dm=mean(dX,3);
[n,k]=size(X);
error.v=(X-repmat(m,n,1));
error.d=dX-repmat(dm,1,1,n);
D=zeros(k*k,size(dX,2));
for i=1:n
    %D=D+kron(error.v(i,:)',eye(k))*error.d(:,:,i)+kron(error.v(i,:)',eye(k));   
    D=D+d_prod(error.v(i,:)',d_trans(error.v(i,:),error.d(:,:,i)),error.v(i,:),error.d(:,:,i));
end    
D=D./(n-1);