function p=d_log_mean_pdf(logp)
M=size(logp.d,1);
p.v=exp(logp.v);
p.d=sparse(1:M, 1:M, p.v)*logp.d;
p.v=mean(p.v);
p.d=mean(p.d);
p.d=1/p.v*p.d;
p.v=log(p.v);