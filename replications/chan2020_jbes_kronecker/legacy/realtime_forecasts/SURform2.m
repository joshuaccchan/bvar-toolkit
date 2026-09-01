%% construct a sparse matrix Xout such that 
% for i=1:r
%    bigX((i-1)*n+1:i*n,:) = kron(speye(n),X(i,:));
% end

function Xout = SURform2( X, n )
repX = kron(X,ones(n,1));
[r c] = size( X );
idi = kron((1:r*n)',ones(c,1));
idj = repmat((1:n*c)',r,1);
Xout = sparse(idi,idj,reshape(repX',n*r*c,1));