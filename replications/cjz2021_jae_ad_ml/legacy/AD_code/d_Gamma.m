function [dG] = d_Gamma(g, alpha) 

  fun = @(t) log(t) .*t.^(alpha-1).*exp(-t);
  num_1 = integral(fun, 0, g)/gamma(alpha);
  
  ps=psi(alpha);
  CDF=gamcdf(g,alpha,1);
  
  if num_1==Inf
  num_1=ps-(1-CDF)*log(g)-meijerG([],[1,1],[0,0,alpha],[],g)/gamma(alpha);
  end
  
  num_2 = ps *CDF;
  
  dG = - (num_1 - num_2) ./ gampdf(g, alpha, 1);
