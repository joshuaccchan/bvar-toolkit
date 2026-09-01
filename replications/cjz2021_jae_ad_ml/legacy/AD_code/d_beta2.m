function [beta,bg, error]=d_beta2(Vg,bg,Z,d)
global X Y BigX K n
 bg.d=d_prod(Vg.v,Vg.d,bg.v,bg.d);
   bg.v=Vg.v*bg.v;

   C=chol(Vg.v,'lower');
   beta.v=bg.v+C*Z;
   beta.d=bg.d+kron(Z',eye(d))*d_cholasky(C,Vg.d);
   
    b=reshape(beta.v,n,K);
   error.v=(Y-X*b'); 
   error.d=-BigX*d_trans(b,beta.d);