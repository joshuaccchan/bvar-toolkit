% This function calculates the variance and covariance for figures 8-9
function [var_4, cov_4] = get_varcov_first2last2(Sig_mat,rev_option)

[T,n,~] = size(Sig_mat);

var_4 = zeros(T,4); cov_4 = zeros(T,4);
for t=1:T
    St = squeeze(Sig_mat(t,:,:));
    var_4(t,:) = [St(1,1),St(2,2),St(19,19),St(20,20)];
    cov_4(t,:) = [St(1,2),St(1,19),St(2,20),St(19,20)];
end
if rev_option == 1
    var_4 = fliplr(var_4); cov_4 = fliplr(cov_4);
end
end