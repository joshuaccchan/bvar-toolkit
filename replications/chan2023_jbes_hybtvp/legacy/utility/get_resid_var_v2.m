% Support function for estimating the hybrid TVP-VAR in Chan (2022)
%
% See:
% Chan, J.C.C. (2022). Large Hybrid Time-Varying Parameter VARs, Journal of
% Business and Economic Statistics, forthcoming

function sig2 = get_resid_var_v2(Y0,Y)
[T,n] = size(Y);
sig2 = zeros(n,1);
tmpY = [Y0(end-4+1:end,:); Y];
for i=1:n
    Z = [ones(T,1) tmpY(4:end-1,:) tmpY(3:end-2,:) tmpY(2:end-3,:) tmpY(1:end-4,:)];
    tmpb = (Z'*Z+1e-4*eye(size(Z,2)))\(Z'*tmpY(5:end,i));
    sig2(i) = mean((tmpY(5:end,i)-Z*tmpb).^2);
end
end