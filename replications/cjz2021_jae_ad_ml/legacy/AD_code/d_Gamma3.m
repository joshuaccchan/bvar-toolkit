% This is written by Dan Zhu(dan.zhu@monash.edu) on the 6th of June 2019.
%This function generates 
%  G ~gamma(alpha,1) 
%and its derivatives

function G=d_Gamma3(alpha,n,k)
G.v=gamrnd(alpha,1,n,k);
dif=(gamcdf(G.v,alpha+1e-8,1)-gamcdf(G.v,alpha-1e-8,1))./2e-8;
%pdf=G.v.^(alpha-1).*exp(-G.v)./gamma(alpha);
pdf=gampdf(G.v,alpha);
G.d=-dif./pdf;