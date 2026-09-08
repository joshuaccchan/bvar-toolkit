% This function obtains the posterior draws of alpha, beta and Sigma
% under the asymmetric conjugate prior
%
% See:
% Chan, J.C.C. (2021). Asymmetric conjugate priors for large Bayesian VARs,
% Quantitative Economics, forthcoming.

function [Alp,Beta,Sig] = sample_ThetaSig(Y0,Y,p,prior,nsim)
[T,n] = size(Y);
tmpY = [Y0(end-p+1:end,:); Y];
Z = zeros(T,n*p); 
for ii=1:p
    Z(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
end
Z = [ones(T,1) Z];
k_beta = n^2*p+n;
k_alp = n*(n-1)/2;
Beta = zeros(nsim,k_beta);
Alp = zeros(nsim,k_alp);
Sig = zeros(nsim,n);
count_alp = 0;
for ii = 1:n
    yi = Y(:,ii);
    ki = n*p+ii;
    mi = [prior.beta0(:,ii);prior.alp0(count_alp+1:count_alp+ii-1)];
    Vi = sparse(1:ki,1:ki,[prior.Vbeta(:,ii);prior.Valp(count_alp+1:count_alp+ii-1)]);
    nui = prior.nu(ii);
    Si = prior.S(ii);    
    Xi = [Z -Y(:,1:ii-1)];
        % compute the parameters of the posterior distribution
    iVi = Vi\speye(ki);
    Kthetai = iVi + Xi'*Xi;
    CKthetai = chol(Kthetai,'lower');    
    thetai_hat = (CKthetai')\(CKthetai\(iVi*mi + Xi'*yi));
    Si_hat = Si + (yi'*yi + mi'*iVi*mi - thetai_hat'*Kthetai*thetai_hat)/2;
        % sample sig and theta
    Sigi = 1./gamrnd(nui+T/2,1./Si_hat,nsim,1);
    U = randn(nsim,ki).*repmat(sqrt(Sigi),1,ki);
    Thetai = repmat(thetai_hat',nsim,1) + U/CKthetai;
    
    Sig(:,ii) = Sigi;
    Beta(:,(ii-1)*(n*p+1)+1:ii*(n*p+1)) =  Thetai(:,1:n*p+1);
    Alp(:,count_alp+1:count_alp+ii-1) =  Thetai(:,n*p+2:end);
    count_alp = count_alp + ii -1;
end   
end