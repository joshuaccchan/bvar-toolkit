function [G,dG]=d_Gamvector2(a,b,da,db,M)


G=gamrnd(a,1,M,1);
for i=1:M
 dG(i,:)=d_Gamma2(G(i,1),a)*da;   
end    
dG=G*db+b*dG;
G=b*G;
dG=dG';