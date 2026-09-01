% This is written by Dan Zhu(dan.zhu@monash.edu) on the 6th of June 2019
% This function computes matrix addition as well as its derivatives.
function C=Maddition(A,B)
if isfield(A,'d')&& isfield(B,'d')==0
    C.v=A.v+B;
    C.d=A.d;
elseif isfield(A,'d')==0&& isfield(B,'d')
    C.v=A+B.v;
    C.d=B.d;
else
C.v=A.v+B.v;
C.d=A.d+B.d;
end