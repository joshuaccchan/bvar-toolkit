% This function computes the forecasts from a large Bayesian VAR with a 
% data-based Minnesota prior 
%
% See:
% Chan, J.C.C. (2021). Minnesota-Type Adaptive Hierarchical Priors for 
% Large Bayesian VARs, International Journal of Forecasting, forthcoming


    % construct a few useful things
tmpyhat1 = zeros(nsim,2*n+1);  %% [point forecasts, prelike, joint den]
tmpyhat4 = zeros(nsim,2*n+1); 
A_id = nonzeros(tril(reshape(1:n^2,n,n),-1)');
A = eye(n);
k_beta = n^2*p+n;
k_alp = n*(n-1)/2;
sig2 = get_resid_var(Y0,Yt);
[C,idx_kappa1,idx_kappa2] = get_C(n,p,sig2);

    % initalize storage
store_kappa = zeros(nsim,2);
store_alp = zeros(nsim,k_alp);
store_beta = zeros(nsim,k_beta);
store_h_T = zeros(nsim,n);  % store only h_T
store_Sigh = zeros(nsim,n);

    % initialize the Markov chain
kappa = [.4,.001,1,100];
h0 = log(sig2);
h = repmat(h0',Tt,1);
Sigh = Sh0;
beta = zeros(k_beta,1);
alp = zeros(k_alp,1);

    % MCMC starts here
for isim = 1:nsim + burnin
        % sample alp and beta    
    [Valp,Vbeta] = getVtheta(idx_kappa1,idx_kappa2,kappa,C,sig2);
    count_alp = 0;
    U = zeros(Tt,n);
    for ii = 1:n
        yi = Yt(:,ii);
        ki = n*p+ii;
        Xi = [Zt -Yt(:,1:ii-1)];        
        iVthetai = sparse(1:ki,1:ki,1./[Vbeta((ii-1)*(n*p+1)+1:ii*(n*p+1));...
            Valp(count_alp+1:count_alp+ii-1)]);
        XiiSighi = Xi'*sparse(1:Tt,1:Tt,exp(-h(:,ii)));
        Kthetai = iVthetai + XiiSighi*Xi;
        CKthetai = chol(Kthetai,'lower');
        thetai_hat = (CKthetai')\(CKthetai\(XiiSighi*yi));
        thetai = thetai_hat + CKthetai'\randn(ki,1);
        U(:,ii) = yi - Xi*thetai;
        
        beta((ii-1)*(n*p+1)+1:ii*(n*p+1)) =  thetai(1:n*p+1);
        alp(count_alp+1:count_alp+ii-1) =  thetai(n*p+2:end);
        count_alp = count_alp + ii-1;
    end
    
        % sample h
    for i=1:n
        Ystar = log(U(:,i).^2 + .0001);
        h(:,i) = SVRW(Ystar,h(:,i),Sigh(i),h0(i));
    end 
    
        % sample kappa1 and kappa2    
    tmpc1 = sum(beta(idx_kappa1).^2./C(idx_kappa1));
    tmpc2 = sum(beta(idx_kappa2).^2./C(idx_kappa2));
    kappa(1) = gigrnd(c01(1)-n*p/2,2*c01(2),tmpc1,1);
    kappa(2) = gigrnd(c02(1)-(n-1)*n*p/2,2*c02(2),tmpc2,1);    
      
        %% sample h0
    Kh0 = sparse(1:n,1:n,1./Sigh + 1./Vh);
    h0_hat = Kh0\(ah./Vh + h(1,:)'./Sigh);
    h0 = h0_hat + chol(Kh0,'lower')'\randn(n,1);
    
        %% sample Sigh
    e = h - [h0';h(1:Tt-1,:)];
    Sigh = 1./gamrnd(nuh0+Tt/2,1./(Sh0 + sum(e.^2)'/2));    
    
    if isim > burnin        
        isave = isim-burnin;        
        store_kappa(isave,:) = kappa(1:2);        
        store_beta(isave,:) =  beta';
        store_alp(isave,:) =  alp';
        store_h_T(isave,:) =  h(end,:)';
        store_Sigh(isave,:) = Sigh';
    end
end
kappa_hat = mean(store_kappa)';
kappaCI = quantile(store_kappa,[0.25 .975]);

    % compute forecasts
for isim = 1:nsim
    alp = store_alp(isim,:)';
    beta = store_beta(isim,:)';
    h_Tp1 = store_h_T(isim,:)';    
    Sigh = store_Sigh(isim,:)';
    
        % trasnform the parameters into reduced-form
    sqrtSigh = sqrt(Sigh);
    A(A_id) = alp;
    h_Tp1 = h_Tp1 + sqrtSigh.*randn(n,1);    
    S = (A\sparse(1:n,1:n,exp(h_Tp1)))/A';
    B = (A\(reshape(beta,n*p+1,n)'))'; 
    xtp1 = [1 reshape(Yt(end:-1:end-p+1,:)',1,n*p)];
    CS = chol(S,'lower');
    for tt=1:4
        EYtp1 = xtp1*B;
        dS = diag(S)';
        if tt == 1
            tmpu = CS\(Y(t+1,:)-EYtp1)';
            lden_joint = -n/2*log(2*pi) -sum(diag(log(CS))) -.5*(tmpu'*tmpu);
            lden = -.5*log(2*pi*dS) - .5*(Y(t+1,:)-EYtp1).^2./dS;
            tmpyhat1(isim,:) = [EYtp1 lden lden_joint];
        elseif tt == 4 && t<=T-tt
            tmpu = CS\(Y(t+4,:)-EYtp1)';
            lden_joint = -n/2*log(2*pi) -sum(diag(log(CS))) -.5*(tmpu'*tmpu);
            lden = -.5*log(2*pi*dS) - .5*(Y(t+4,:)-EYtp1).^2./dS;
            tmpyhat4(isim,:) = [EYtp1 lden lden_joint];
        end
        Ytp1 = EYtp1 + (CS*randn(n,1))';
        xtp1 = [1 Ytp1 xtp1(2:end-n)];
        
        h_Tp1 = h_Tp1 + sqrtSigh.*randn(n,1);    
        S = (A\sparse(1:n,1:n,exp(h_Tp1)))/A';
        CS = chol(S,'lower');
    end
    
end


