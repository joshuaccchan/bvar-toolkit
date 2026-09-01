function dC=d_dotinv(A,dA)
Da=-diag(A.^2);
dC=Da*dA;
