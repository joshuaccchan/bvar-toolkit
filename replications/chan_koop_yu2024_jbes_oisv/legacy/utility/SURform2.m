function X_out = SURform2(X, n)
    repX = kron(X,ones(n,1));
    [r,c] = size( X );
    idi = kron((1:r*n)',ones(c,1));
    idj = repmat((1:n*c)',r,1);
    X_out = sparse(idi,idj,reshape(repX',n*r*c,1));
end