function D=d_trace(dA,n)
ind=1+(n+1)*(0:n-1);
D=sum(dA(ind,:));