% This function constructs the time-varying covariance matrices 
function Sigt = construct_Sigt(h,B0)
[T,n] = size(h);
Sigt = zeros(T,n,n);
for t=1:T
    Sigt(t,:,:) = (B0\sparse(1:n,1:n,exp(h(t,:))))/(B0');
end
end