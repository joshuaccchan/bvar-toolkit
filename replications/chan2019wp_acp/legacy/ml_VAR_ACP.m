% This script computes the log marginal likelihood of the BVAR with the
% asymmetric conjugate prior
%
% See:
% Chan, J.C.C. (2019). Asymmetric conjugate priors for large Bayesian VARs,
% CAMA Working Papers 51/2019

function lml = ml_VAR_ACP(p,Y,Y0,Z,kappa)
    [T,n] = size(Y);
    sig2 = get_resid_var(Y0,Y);
    lml = -n*T/2*log(2*pi);
    for i = 1:n        
        yi = Y(:,i);
        ki = n*p+i;
        [mi,Vi,nui,Si] = prior_ACPi(n,p,i,kappa,sig2);
        Xi = [Z -Y(:,1:i-1)];
        
        iVi = Vi\speye(ki);
        Kthetai = iVi + Xi'*Xi;
        CKthetai = chol(Kthetai,'lower');    
        thetai_hat = CKthetai'\(CKthetai\(iVi*mi + Xi'*yi));        
        Si_hat = Si + (yi'*yi + mi'*iVi*mi - thetai_hat'*Kthetai*thetai_hat)/2;
        
        lml = lml -1/2*(sum(log(diag(Vi))) + 2*sum(log(diag(CKthetai)))) ...
            + nui*log(Si) - (nui+T/2)*log(Si_hat) + gammaln(nui+T/2) - gammaln(nui);        
    end    
end