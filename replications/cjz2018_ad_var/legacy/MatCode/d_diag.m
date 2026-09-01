function [ddiagX] = d_diag(dX, k)  
s=k*(0:k-1)+(1:k);
In=sparse(s,1:k,ones(k,1));
ddiagX =In * dX;

