function MLCE= SF_CE_AD(store_a,store_Sig,store_Omega,Y,prior,R,Ind,Ind_a)
R = 20*ceil(R/20);   % make R divisible by 20
na = size(store_a.v,2);
n = size(Y,2); r = size(store_Omega.v,2);

    % obtain parameters for the IS density
[a_bar,Cov,CDa_bar,...
    nusig2_bar,Ssig2_bar,...
    nuomega2_bar,Somega2_bar]=IS_parameters(store_a,store_Sig,store_Omega, R);

    % obtain IS draws from the optimal density
[invSig_IS,invOmega_IS,mA]=Sample_IS(a_bar,CDa_bar,...
    nusig2_bar,Ssig2_bar,...
    nuomega2_bar,Somega2_bar,R);

Sig_IS=matrix_dot_division(ones(n,R),invSig_IS);
Omega_IS=matrix_dot_division(ones(r,R),invOmega_IS);
% obtain the prior
Ssig2_bar=matrix_dot_division(ones(n,1),Ssig2_bar);
Somega2_bar=matrix_dot_division(ones(r,1),Somega2_bar);
prior=prior(mA,Sig_IS,Omega_IS);
    % obtain the IS density
g_IS =G_IS(mA,Sig_IS,Omega_IS,...
    a_bar,CDa_bar,nusig2_bar,Ssig2_bar,nuomega2_bar,Somega2_bar);

store_w.v = zeros(R,1);
store_w.d=zeros(R,4);
K_rn=commutation_matrix(r,n);
K_nn=commutation_matrix(n,n);
E_n=elimination_matrix(n);
for isim  = 1:R
    invSig.v = invSig_IS.v(:,isim);
    invSig.d=invSig_IS.d((isim-1)*n+1:isim*n,:); 
    invOmega.v = invOmega_IS.v(:,isim); 
    invOmega.d=invOmega_IS.d((isim-1)*r+1:isim*r,:);  
    a.v=mA.v(isim,:);
    a.d=mA.d(isim+(0:na-1)*R,:);
    llike = Intlike_AD(Y,a,invSig,invOmega,Ind,Ind_a,K_rn,K_nn,E_n);  
    store_w.v(isim)=llike.v;
    store_w.d(isim,:)=llike.d;
end
MLCE=d_CE(prior,store_w,g_IS);
end