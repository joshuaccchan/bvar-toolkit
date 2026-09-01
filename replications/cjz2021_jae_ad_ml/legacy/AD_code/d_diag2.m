% This is written by Dan Zhu(dan.zhu@monash.edu) on the 16th of June 2018
% It computes derivative of the operation 'diag'. When vector==true, it
% turn vector to diagonal matrix. 

function C=d_diag2(X,vector,k)
s=k*(0:k-1)+(1:k);
In=sparse(s,1:k,ones(k,1));
if vector==true
C.d =sparse(In * X.d);
C.v=sparse(1:k,1:k,X.v);
else   
C.d =In' * X.d;
C.v=diag(X.v);
end    
    
    
    