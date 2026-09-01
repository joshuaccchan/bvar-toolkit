function [Omega,invOmega, S2]=Sample_Omega(F,Somega2,nuomega2,T,r)
d=size(F.d,2);
 E2=Sum_matrix(matrix_dot_times(F,F),2);
    S2.v=0.5*E2.v+Somega2; S2.d=0.5*E2.d;S2.d(:,4)=S2.d(:,4)+1;
    Omega=d_Gamma3(nuomega2+T/2,r,1);
    Omega.d=[sparse(r,3),Omega.d,sparse(r,d-4)];
    Omega =matrix_dot_division(S2,Omega);
    invOmega=matrix_dot_division(ones(r,1),Omega);