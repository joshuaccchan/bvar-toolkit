function d=d_Incomplete_Gamma(alpha,x)

d=log(x)*gamcdf(x,alpha,1)*gamma(alpha)+x*meijerG(3,[],[],alpha,x)