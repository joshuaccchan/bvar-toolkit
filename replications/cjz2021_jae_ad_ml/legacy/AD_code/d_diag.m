%This function takes in derivatives of a vector to get the derivatives of
%its corresponding diagonal matrix

function [ddiagX] = d_diag(dX, k)  
s=k*(0:k-1)+(1:k);
In=sparse(s,1:k,ones(k,1));
ddiagX =sparse(In * dX);

