function [F]=Sample_F(invSig,invOmega,A,K_rn,K_rr,E_r,I_r,I_rr,r,T, Y)
 AiSig = Mtimes(A,invSig);
 At=Mtrans(A,K_rn);
    AiSigA = Mtimes(AiSig,At);
  
    K=Maddition(invOmega,AiSigA);
    mu=Mtimes(M_division(AiSig,K,true),Y');
    L=Cholasky(K,K_rr,E_r,I_r,I_rr);
    Lt=Mtrans(L,K_rr);
    F=Maddition(mu,M_division(randn(r,T),Lt,true));