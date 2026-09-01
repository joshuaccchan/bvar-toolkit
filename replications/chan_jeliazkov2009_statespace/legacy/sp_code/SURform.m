% % =======================================================================
% % This support function takes a T x k matrix and returns a T x Tk 
% %     sparse matrix in the form of a seemingly unrelated regression
% %
% % Please report any errors to joshuacc.chan@gmail.com
% % =======================================================================
function Xout = SURform(X)
[r c] = size( X );
idi = kron((1:r)',ones(c,1));
idj = (1:r*c)';
Xout = sparse(idi,idj,reshape(X',r*c,1));