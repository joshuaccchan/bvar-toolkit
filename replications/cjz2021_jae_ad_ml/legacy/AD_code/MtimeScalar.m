% This is written by Dan Zhu(dan.zhu@monash.edu) on the 6th of June 2019
% This function computes matrix times a scalor as well as its derivatives.
function C=MtimeScalar(A,c)
if isfield(A,'d')&& isfield(c,'d')==0
    C.v=A.v.*c;
    C.d=A.d*c;
elseif isfield(A,'d')==0&& isfield(c,'d')
    C.v=A*c.v;
    C.d=A(:)*c.d;
else
C.v=A.v*c.v;
C.d=A.v(:)*c.d+ c.v*A.d;
end

