function [dG] = d_Gamma2(g, alpha) 

 num1=(gamcdf(g,alpha+1e-8,1)-gamcdf(g,alpha-1e-8,1))./2e-8;
  
  dG = - num1./ gampdf(g, alpha, 1);
