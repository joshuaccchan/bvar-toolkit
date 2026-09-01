function dC=d_dotproduct(A,dA, B,dB)

Da=diag(A);
Db=diag(B);

      dC =Da*dB +Db*dA;
