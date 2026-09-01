% Support function for estimating the hybrid TVP-VAR in Chan (2022)
%
% See:
% Chan, J.C.C. (2022). Large Hybrid Time-Varying Parameter VARs, Journal of
% Business and Economic Statistics, forthcoming

function [gami,thetai,Ui,tilde_thetai,lp_gami] = sample_gam_thetai_ver2(Xi,thetai0,...
    Sigthetai,Yi,hi,p0i,k,is_gamfixed,gami)
[T,ki] = size(Xi);
iehi = exp(-hi);
Xithetai0 = Xi*thetai0;
Ei = Yi-Xithetai0;
m = ki-k;
lp_gami = zeros(4,1);
if ~is_gamfixed % if gam^theta_i is not fixed, sample gam_thetai marginal of thetai
    pri_gam = [1-p0i(1);1-p0i(1);p0i(1);p0i(1)].*[1-p0i(2);p0i(2);1-p0i(2);p0i(2)];
    c1 = -T/2*log(2*pi) -.5*sum(hi);
    c2 = -.5*iehi'*(Ei.^2);
        % likelihood value for gam = (0,0)
    llike00 = c1 + c2;
    
        % compute likelihood value for gam = (0,1)
    if m > 0
        H_alpi = speye(T*m) - sparse(m+1:T*m,1:(T-1)*m,ones((T-1)*m,1),T*m,T*m);
        Z01 = SURform(Xi(:,k+1:end).*repmat(sqrt(Sigthetai(k+1:end))',T,1));
        Z01Sig = Z01'*sparse(1:T,1:T,iehi);
        HH_alpi = H_alpi'*H_alpi;
        Kalpi = HH_alpi + Z01Sig*Z01;
        chol_Kalpi = chol(Kalpi,'lower');
        alpi_hat = chol_Kalpi'\(chol_Kalpi\(Z01Sig*Ei));
        llike01 = c1 - sum(log(diag(chol_Kalpi))) + c2 + .5*alpi_hat'*Kalpi*alpi_hat;
    else
        llike01 = min(-1e10,llike00-100);
    end
    
        % compute likelihood value for gam = (1,0)
    H_betai = speye(T*k) - sparse(k+1:T*k,1:(T-1)*k,ones((T-1)*k,1),T*k,T*k);
    Z10 = SURform(Xi(:,1:k).*repmat(sqrt(Sigthetai(1:k))',T,1));
    Z10Sig = Z10'*sparse(1:T,1:T,iehi);
    HH_betai = H_betai'*H_betai;
    Kbetai = HH_betai + Z10Sig*Z10;
    chol_Kbetai = chol(Kbetai,'lower');
    betai_hat = chol_Kbetai'\(chol_Kbetai\(Z10Sig*Ei));
    llike10 = c1 - sum(log(diag(chol_Kbetai))) + c2 + .5*betai_hat'*Kbetai*betai_hat;
        
        % compute likelihood value for gam = (1,1) 
    if m > 0
        H_thetai = speye(T*ki) - sparse(ki+1:T*ki,1:(T-1)*ki,ones((T-1)*ki,1),T*ki,T*ki);
        Z11 = SURform(Xi.*repmat(sqrt(Sigthetai)',T,1));
        Z11Sig = Z11'*sparse(1:T,1:T,iehi);
        HH_thetai = H_thetai'*H_thetai;
        Kthetai = HH_thetai + Z11Sig*Z11;
        chol_Kthetai = chol(Kthetai,'lower');
        thetai_hat = chol_Kthetai'\(chol_Kthetai\(Z11Sig*Ei));
        llike11 = c1 - sum(log(diag(chol_Kthetai))) + c2 + .5*thetai_hat'*Kthetai*thetai_hat;
    else
        llike11 = min(-1e10,llike10-100);
    end
    
    llike = [llike00,llike01,llike10,llike11]';
    maxllike = max(llike);
    prob_gami_ker = exp(llike-maxllike).*pri_gam;
    prob_gami = prob_gami_ker/sum(prob_gami_ker);    
    lp_gami = llike-maxllike+log(pri_gam)-log(sum(prob_gami_ker));    
    tmpdraw = find(cumsum(prob_gami) > rand,1);
    if tmpdraw == 1
        gami = [0,0];
    elseif tmpdraw == 2
        gami = [0,1];
    elseif tmpdraw == 3
        gami = [1,0];
    elseif tmpdraw == 4
        gami = [1,1];
    end        
else % if gam^theta_i is fixed, compute a few things for sampling thetai
    if gami(1) == 1 && gami(2) == 1
        H_thetai = speye(T*ki) - sparse(ki+1:T*ki,1:(T-1)*ki,ones((T-1)*ki,1),T*ki,T*ki);
        Z11 = SURform(Xi.*repmat(sqrt(Sigthetai)',T,1));
        Z11Sig = Z11'*sparse(1:T,1:T,iehi);
        HH_thetai = H_thetai'*H_thetai;
        Kthetai = HH_thetai + Z11Sig*Z11;
        chol_Kthetai = chol(Kthetai,'lower');
        thetai_hat = chol_Kthetai'\(chol_Kthetai\(Z11Sig*Ei)); 
    elseif gami(1) == 1 &&  gami(2) == 0
        H_betai = speye(T*k) - sparse(k+1:T*k,1:(T-1)*k,ones((T-1)*k,1),T*k,T*k);
        Z10 = SURform(Xi(:,1:k).*repmat(sqrt(Sigthetai(1:k))',T,1));
        Z10Sig = Z10'*sparse(1:T,1:T,iehi);
        HH_betai = H_betai'*H_betai;
        Kbetai = HH_betai + Z10Sig*Z10;
        chol_Kbetai = chol(Kbetai,'lower');
        betai_hat = chol_Kbetai'\(chol_Kbetai\(Z10Sig*Ei));
    elseif gami(1) == 0 && gami(2) == 1
        H_alpi = speye(T*m) - sparse(m+1:T*m,1:(T-1)*m,ones((T-1)*m,1),T*m,T*m);
        Z01 = SURform(Xi(:,k+1:end).*repmat(sqrt(Sigthetai(k+1:end))',T,1));
        Z01Sig = Z01'*sparse(1:T,1:T,iehi);
        HH_alpi = H_alpi'*H_alpi;
        Kalpi = HH_alpi + Z01Sig*Z01;
        chol_Kalpi = chol(Kalpi,'lower');
        alpi_hat = chol_Kalpi'\(chol_Kalpi\(Z01Sig*Ei));   
    end
end

    % sample thetai
mu_thetai = kron(ones(T,1),thetai0);
if gami(1) == 1 &&  gami(2) == 1
    tilde_thetai = thetai_hat + chol_Kthetai'\randn(T*ki,1);
    thetai = mu_thetai + repmat(sqrt(Sigthetai),T,1).*tilde_thetai;    
    Ui = Ei - Z11*tilde_thetai;
elseif gami(1) == 1 &&  gami(2) == 0   
    tilde_betai = betai_hat + chol_Kbetai'\randn(T*k,1);
    tilde_Betai = reshape(tilde_betai,k,T)';
    tilde_thetai = reshape([tilde_Betai sparse(T,m)]',T*ki,1);  % tilde_alpi isn't used so it's set to 0
    thetai = mu_thetai + repmat(sqrt(Sigthetai),T,1).*tilde_thetai;     
    Ui = Ei - Z10*tilde_betai;    
elseif gami(1) == 0 &&  gami(2) == 1
    tilde_alpi = alpi_hat + chol_Kalpi'\randn(T*m,1);
    tilde_Alpi = reshape(tilde_alpi,m,T)';
    tilde_thetai = reshape([sparse(T,k), tilde_Alpi]',T*ki,1);
    thetai = mu_thetai + repmat(sqrt(Sigthetai),T,1).*tilde_thetai; 
    Ui = Ei - Z01*tilde_alpi; 
elseif gami(1) == 0 &&  gami(2) == 0
    thetai = mu_thetai;
    Ui = Ei;    
    tilde_thetai = sparse(T*ki,1); % tilde_thetai isn't used so it's set to 0
end

end