% This script estimates the standard BVAR - no MCMC is needed
%
% See:
% Chan, J.C.C. (2020). Large Bayesian VARs: A flexible Kronecker error 
% covariance structure, Journal of Business and Economic Statistics, 
% 38(1), 68-79.

%% prior
S0 = eye(n);
nu0 = n+3;
construct_prior_A;

%% compute the parameters for the posterior density
X = zeros(T,n*p); 
for i=1:p
    X(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end
X = [ones(T,1) X];
XX = X'*X;
Atilde = XX\(X'*shortY);
KA = sparse(1:k,1:k,1./VA0) + XX;
    % posterior mean of the VAR coefficients, arranged as a k by n matrix
Ahat = KA\(sparse(1:k,1:k,VA0)\A0 + XX*Atilde); 
Shat = S0 + A0'*sparse(1:k,1:k,1./VA0)*A0 + shortY'*shortY ...
    - Ahat'*KA*Ahat;
Shat = (Shat+Shat')/2;

figure;
colormap('hsv');
imagesc(Ahat);
colorbar;
box off;
title('Heat map of the VAR coefficients');    

%% compute the marginal likelihood if cp_ml = 1
if cp_ml
    A = Ahat;
    Sig = Shat/(T+nu0);
    CSig = chol(Sig,'lower');
    tmp = (shortY-X*A)/CSig';
    s2 = sum(tmp.^2,2);
    llike =  -T*n/2*log(2*pi) - T*sum(log(diag(CSig))) -.5*sum(s2);
    ML = llike + lniwpdf(A,Sig,A0,sparse(1:k,1:k,1./VA0),nu0,S0) ...
        - lniwpdf(A,Sig,Ahat,KA,nu0+T,Shat);
    
    fprintf('log marginal likelihood of BVAR: %.1f \n', ML); 
end

