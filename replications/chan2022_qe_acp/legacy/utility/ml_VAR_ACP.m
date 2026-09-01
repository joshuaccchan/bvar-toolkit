% This function computes the marginal likelihood value under the asymmetric
% conjugate prior
%
% See:
% Chan, J.C.C. (2022). Asymmetric Conjugate Priors for Large Bayesian VARs,
% Quantitative Economics, 13(3): 1145-1169

function lml = ml_VAR_ACP(p,Y,Z,prior)
    [T,n] = size(Y);    
    lml = -n*T/2*log(2*pi);
    count_alp = 0;
    for ii = 1:n        
        yi = Y(:,ii);
        ki = n*p+ii;
        mi = [prior.beta0(:,ii);prior.alp0(count_alp+1:count_alp+ii-1)];
        Vi = sparse(1:ki,1:ki,[prior.Vbeta(:,ii);prior.Valp(count_alp+1:count_alp+ii-1)]);
        nui = prior.nu(ii);
        Si = prior.S(ii);    
        Xi = [Z -Y(:,1:ii-1)];
        
        iVi = Vi\speye(ki);
        Kthetai = iVi + Xi'*Xi;
        CKthetai = chol(Kthetai,'lower');    
        thetai_hat = CKthetai'\(CKthetai\(iVi*mi + Xi'*yi));        
        Si_hat = Si + (yi'*yi + mi'*iVi*mi - thetai_hat'*Kthetai*thetai_hat)/2;
        
        lml = lml -1/2*(sum(log(diag(Vi))) + 2*sum(log(diag(CKthetai)))) ...
            + nui*log(Si) - (nui+T/2)*log(Si_hat) + gammaln(nui+T/2) - gammaln(nui); 
        
        count_alp = count_alp + ii-1;        
    end    
end