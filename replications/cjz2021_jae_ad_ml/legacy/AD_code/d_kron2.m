% This is written by Dan Zhu(dan.zhu@monash.edu) on the 18th of Jan 2020
% for computing the kronecker product of two matrices and the assocaited
% derivatives
function C = d_kron(A, B,K_qm)
if isfield(A,'d')&&isfield(B,'d')==0
   [m,n] = size(A.v);
  [p,q] = size(B); 
  if isempty(K_qm)
  K_qm = commutation_matrix(q, m);
  end
  f1 = sparse(kron(kron(speye(n), K_qm), speye(p)));
  C.v=sparse(kron(A.v,B));
  C.d = f1 * sparse(kron(A.d, B(:)));
elseif isfield(A,'d')==0&&isfield(B,'d')
      [m,n] = size(A);
  [p,q] = size(B.v);
   if isempty(K_qm)
  K_qm = commutation_matrix(q, m);
   end
  f1 = sparse(kron(kron(speye(n), K_qm), speye(p)));
  C.v=sparse(kron(A,B.v));
      C.d = f1 * sparse(kron(A(:), B.d));
else  
  [m,n] = size(A.v);
  [p,q] = size(B.v);
  if isempty(K_qm)
  K_qm = commutation_matrix(q, m);
  end
  f1 = sparse(kron(kron(speye(n), K_qm), speye(p)));
  C.v=sparse(kron(A.v,B.v));
      C.d = f1 * (sparse(kron(A.v(:), B.d)) + sparse(kron(A.d, B.v(:))));
end




