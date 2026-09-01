function dC=d_dotinv(invA,dA)
Da=-diag(invA.^2);
dC=Da*dA;
