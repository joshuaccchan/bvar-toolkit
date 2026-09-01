
function S=d_gamMLE(a,b,S,dS,n)      
   
V1=1./S; V2=1./S.^2;
H=[-n*psi(1,a),-n/b;
   -n/b,     a*n/b^2-2/b^3*sum(V1)];
 
       Delta_s=[-V1;
           -1/b^2*V2];
       S=-H\Delta_s*dS;
end