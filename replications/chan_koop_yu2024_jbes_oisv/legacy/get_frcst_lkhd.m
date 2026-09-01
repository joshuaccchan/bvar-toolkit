% This function calculates the rmsfe and likelihood
function [rmsfe, alpl] = get_frcst_lkhd(yhat,option)

n = 20;
for i=1:12
    yhati = yhat(:,:,i);
    rmsfe(:,i)=sqrt(mean((yhati(13-i:end-(i-1),1:n)-yhati(13-i:end-(i-1),n+1:2*n)).^2));
    alpl(:,i)=mean(yhati(13-i:end-(i-1),2*n+1:3*n));
end

if option == 1 % order 2
    rmsfe=flipud(rmsfe);
    alpl=flipud(alpl);
end

end