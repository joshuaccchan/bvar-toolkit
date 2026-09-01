function [dG] = d_Gamma(g, alpha) 
  fun = @(t) log(t) .* gampdf(t, alpha, 1);
  num_1 = integral(fun, 0, g);
  num_2 = psi(alpha) * gamcdf(g, alpha, 1);
  dG = - (num_1 - num_2) ./ gampdf(g, alpha, 1);
