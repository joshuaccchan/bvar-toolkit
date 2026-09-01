function [Sig,invSig,S]=Sample_Sig(A,Ft,Y,Ssig2,T,nusig2,n)
  E=Msubtraction(Y,Mtimes(Ft,A));
    E=Sum_matrix(matrix_dot_times(E,E),1);
    S.v=0.5*E.v'+Ssig2;
    S.d=0.5*E.d;S.d(:,3)=S.d(:,3)+1;
    Sig=d_Gamma3(nusig2+T/2,n,1);
    Sig.d=[sparse(n,2),Sig.d,sparse(n,size(S.d,2)-3)];
    Sig=matrix_dot_division(S,Sig);
    invSig=matrix_dot_division(ones(n,1),Sig);