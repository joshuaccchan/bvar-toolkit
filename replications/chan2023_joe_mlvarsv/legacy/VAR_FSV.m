% This script estimates the VAR-FSV model and computes the associated 
% marginal likelihood
%
% See:
% Chan, J.C.C. (2022). Comparing Stochastic Volatility Specifications for 
% Large Bayesian VARs, Journal of Econometrics, forthcoming.

    % initialize storage
store_l = zeros(nsim,n*r-r*(r+1)/2);
store_A = zeros(k,n);
store_h = zeros(nsim,T,n+r);
store_F = zeros(nsim,T,r);
store_hpara = zeros(nsim,3*(n+r)); % [mu' phi sig2];
count_phi = zeros(n+r,1);
store_kappa = zeros(nsim,2);
L_idx = find(tril(ones(n,r),-1)~=0); % index of free elements of L

    % initialize the Markov chain
A = (X'*X + sparse(1:k,1:k,1./[100,.1*ones(1,k-1)]))\(X'*Y);
L = [eye(r); ones(n-r,r)];
E = Y-X*A;
varY = var(Y)'; 
mu = [log(varY/2); log(mean(varY))*ones(r,1)];
phi = .9 + .09*rand(n+r,1);
sig2 = .05 + .05*rand(n+r,1);
Uy = Y-X*A;
h = [zeros(T,n) repmat(mu(n+1:end)',T,1)];
for ii=1:n
    h(:,ii) = getARh_approx1N(Uy(:,ii).^2,mu(ii),phi(ii),sig2(ii));
end
[C_alp,idx_kappa1,idx_kappa2] = get_C(n,p,sig2_hat);

    % MCMC starts here
randn('seed',sum(clock*100)); rand('seed',sum(clock*1000));
start_time = clock;
for isim = 1:nsim + burnin    
        % sample F    
    e = reshape((Y-X*A)',T*n,1);
    Xf = kron(speye(T),L);
    XfiSig = Xf'*sparse(1:T*n,1:T*n,reshape(exp(-h(:,1:n))',T*n,1)); 
    Kf = sparse(1:T*r,1:T*r,reshape(exp(-h(:,n+1:end))',T*r,1)) + XfiSig*Xf;
    f_hat = Kf\(XfiSig*e); 
    f = f_hat + chol(Kf,'lower')'\randn(T*r,1);
    F = reshape(f,r,T)';
    
        % sample L and A
    [~,Hyper.Valp] = prior_Minn(p,kappa1,kappa2,kappa3,Y0,Y);
    for ii = 1:n 
        if ii <= r
            Zi = [X F(:,1:ii-1)];
            iVthetai = sparse(1:k+ii-1,1:k+ii-1,1./[Hyper.Valp((ii-1)*k+1:ii*k); ...
                Hyper.Vl*ones(ii-1,1)]);
            thetai0 = [Hyper.alp0((ii-1)*k+1:ii*k); Hyper.l0*ones(ii-1,1)];
            ZiSigi = Zi'*sparse(1:T,1:T,exp(-h(:,ii))); 
            dthetai = iVthetai*thetai0 + ZiSigi*(Y(:,ii)-F(:,ii)); 
        else
            Zi = [X F];            
            iVthetai = sparse(1:k+r,1:k+r,1./[Hyper.Valp((ii-1)*k+1:ii*k); ...
                Hyper.Vl*ones(r,1)]);
            thetai0 = [Hyper.alp0((ii-1)*k+1:ii*k); Hyper.l0*ones(r,1)];
            ZiSigi = Zi'*sparse(1:T,1:T,exp(-h(:,ii)));
            dthetai = iVthetai*thetai0 + ZiSigi*Y(:,ii);
        end
        Kthetai = iVthetai + ZiSigi*Zi;
        CKthetai = chol(Kthetai,'lower');
        thetai_hat = CKthetai'\(CKthetai\dthetai);            
        thetai = thetai_hat + CKthetai'\randn(k+min(ii-1,r),1);                    
        
        A(:,ii) = thetai(1:k);
        L(ii,1:min(ii-1,r)) = thetai(k+1:end);
    end
    alp = A(:);
          
        % sample h        
    Uy = Y -X*A -F*L';    
    ystar = log([Uy F].^2 + 1e-4);
    for jj=1:n+r
        h(:,jj) = sample_SV(ystar(:,jj),h(:,jj),mu(jj),phi(jj),sig2(jj));
    end     
   
        % sample mu, phi and sig2
    [mu,phi,sig2,flag_phi] = sample_SVpara(h,mu,phi,Hyper);
    count_phi = count_phi + flag_phi;   
    
        % sample kappa1 and kappa2
    if is_kappafixed
        % do nothing
    elseif is_kappasym              
        tmpc1 = sum(alp(idx_kappa1).^2./C_alp(idx_kappa1)) + ...
            sum(alp(idx_kappa2).^2./C_alp(idx_kappa2));
        kappa1 = gigrnd(Hyper.c0(1,1)-n^2*p/2,2*Hyper.c0(1,2),tmpc1,1);
        kappa2 = kappa1;        
    else
        tmpc1 = sum(alp(idx_kappa1).^2./C_alp(idx_kappa1));
        tmpc2 = sum(alp(idx_kappa2).^2./C_alp(idx_kappa2));    
        kappa1 = gigrnd(Hyper.c0(1,1)-n*p/2,2*Hyper.c0(1,2),tmpc1,1);
        kappa2 = gigrnd(Hyper.c0(2,1)-(n-1)*n*p/2,2*Hyper.c0(2,2),tmpc2,1);    
    end
    
    if isim > burnin
        isave = isim - burnin; 
        store_l(isave,:) = L(L_idx);
        store_A = store_A + A;
        store_h(isave,:,:) = h;        
        store_F(isave,:,:) = F;
        store_hpara(isave,:) =  [mu',phi',sig2'];        
        store_kappa(isave,:) = [kappa1,kappa2];
    end
    
    if ( mod(isim, 10000) ==0 )
        disp(  [ num2str(isim) ' loops... ' ] )
    end    
end
disp( ['Estimating ' model_name ' takes '  num2str( etime( clock, start_time) ) ' seconds' ] );

l_mean = mean(store_l);
L_mean = [eye(r);ones(n-r,r)]; L_mean(L_idx) = l_mean;
A_mean = store_A/nsim;
h_mean = squeeze(mean(store_h));
F_mean = squeeze(mean(store_F));
hpara_mean = mean(store_hpara)';
kappa_mean = mean(store_kappa)';

if cp_ml
    disp(['Computing the marginal likelihood of ' model_name '... ']);
    start_time = clock;    
    [lml,lmlstd,store_w] = ml_var_fsv(X,Y,Y0,M,Hyper,flag_marg,store_h,...
        store_hpara,store_l,store_kappa,is_kappafixed,is_kappasym);
    disp( ['ML computation takes '  num2str(etime(clock, start_time)) ' seconds']);
end
disp(' ' );


