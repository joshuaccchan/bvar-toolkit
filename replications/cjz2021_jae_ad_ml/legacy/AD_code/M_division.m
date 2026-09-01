% This is written by Dan Zhu(dan.zhu@monash.edu) on the 6th of Sept 2019 
%for computing derivatives
% of matrix division, i.e. if pre==true(B\A), else  A/B;

function C=M_division(A,B,pre)
if pre==true
if isfield(A,'d')&&isfield(B,'d')==0
    invB=speye(size(B,2))/B;
    C.d=kron(speye(size(A.v,2)),invB) * A.d;
    C.v=invB*A.v;
elseif isfield(A,'d')==0&&isfield(B,'d')
    C.v=B.v\A;
    C.d=-sparse(kron(C.v',speye(size(B.v,2))/B.v))*B.d;   
else  
 invB=speye(size(B.v,2))/B.v;
 C.v=invB*A.v;
 C.d=kron(speye(size(A.v,2)),invB) * (A.d-sparse(kron(C.v',speye(size(B.v,2))))*B.d);
end
else
   
    if isfield(A,'d')&&isfield(B,'d')==0
    invB=speye(size(B,2))/B';
    C.d=sparse(kron(invB,speye(size(A.v,1)))) * A.d;
    C.v=A.v*invB';
elseif isfield(A,'d')==0&&isfield(B,'d')
    C.v=A/B.v;
    C.d=-sparse(kron(speye(size(B.v,2))/B.v'),C.v)*B.d;   
    else  

 invB=speye(size(B.v,2))/B.v';
 C.v=A.v*invB';
 C.d=kron(invB,speye(size(A.v,1)))*(A.d-sparse(kron(speye(size(B.v,2)),C.v))*B.d);

    end
end   
end

