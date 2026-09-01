

function lpa=logmean(store_lpa,nsim)
lpa.v=exp(store_lpa.v);
lpa.d=sparse(1:nsim,1:nsim,lpa.v)*store_lpa.d;
lpa.v=mean(lpa.v);
lpa.d=mean(lpa.d);
lpa.d=lpa.d./lpa.v;lpa.v=log(lpa.v);