% this function takes in derivatives of a matrix to the derivatives of its diagonal terms 

function [ddiagX] = d_diagg(dX, k) 
s=k*(0:k-1)+(1:k);
In=sparse(s,1:k,ones(k,1));
ddiagX =In' * dX;