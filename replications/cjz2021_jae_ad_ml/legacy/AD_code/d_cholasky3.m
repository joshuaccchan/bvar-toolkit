function d=d_cholasky3(U,dC,n,k, upper)
d=sparse(n^2,k);
inv=U\speye(n);
invt=inv';
if upper==false
for i=1:k
   C=inv*reshape(dC(:,i),n,n)*invt;
   C=U*(tril(C,1)+0.5*diag(diag(C)));
   d(:,i)=C(:);
end    
else
   for i=1:k
   C=sparse(invt*reshape(dC(:,i),n,n)*inv);
   C=(triu(C,1)+0.5*diag(diag(C)))*U;
   d(:,i)=C(:);
   end
end    
    
