% This is written by Dan Zhu(dan.zhu@monash.edu) on the 16th of Jan 2020
% This function computes matrix substraction as well as its derivatives.

function C=Msubtraction(A,B)

if isfield(A,'d')&&isfield(B,'d')==0
    C.v=A.v-B;
    C.d=A.d;
    
elseif isfield(A,'d')==0&&isfield(B,'d')
    C.v=A-B.v;
    C.d=-B.d;
else
C.v=A.v-B.v;
C.d=A.d-B.d;
end