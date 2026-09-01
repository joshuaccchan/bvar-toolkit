% This is written by Dan Zhu(dan.zhu@monash.edu) on the 6th of June 2019
% This function computes matrix inverse as well as its derivatives.
function [inv_A] = Minverse(A)
dim=size(A.v,1);
inv_A.v=A.v\speye(dim);
f1 = sparse(kron(inv_A.v', inv_A.v));
 inv_A.d = -f1 * A.d;
