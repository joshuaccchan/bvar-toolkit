var_small = [1 8 12 14]'; % real GDP growth, unemployment, PCE inflation, Fed Funds rate
n = length(var_small);
k_beta = n*(n-1)/2;      % dimension of the impact matrix
k_alp = n^2*p + n;       % dimension of states 
Y0 = data_t(1:4,var_small);
shortYt = data_t(5:end,var_small);
yt = reshape(shortYt',Tt*n,1);

%% prior
ah = zeros(n,1); Vh = 10*ones(n,1);
nuh0 = 5*ones(n,1); Sh0 = .01*ones(n,1).*(nuh0-1);
beta0 = zeros(k_beta,1); 
alp0 = zeros(k_alp,1); 
c1 = .2^2; c2 = 100; c3 = 10;
Valp = zeros(n*p+1,1);
sig2 = zeros(n,1);
    % construct a prior
tmpY = [Y0(end-p+1:end,:); shortYt];
for i=1:n
    Z = [ones(Tt,1) tmpY(4:end-1,i) tmpY(3:end-2,i) tmpY(2:end-3,i)...
        tmpY(1:end-4,i)];
    tmpb = (Z'*Z)\(Z'*tmpY(5:end,i));
    sig2(i) = mean((tmpY(5:end,i)-Z*tmpb).^2);
end
for i=1:n*p+1
    l = ceil((i-1)/n); 
    idx = mod(i-1,n); % variable index
    if idx==0
        idx = n;
    end        
    if i==1 % intercept
        Valp(1) = c2;     
    else
        Valp(i) = c1/(l^2*sig2(idx));        
    end   
end
Valp = repmat(Valp,n,1);
Vbeta = c3*ones(k_beta,1);

%% compute and define a few things
tmpY = [Y0; shortYt];
x = zeros(Tt,n*p); 
for i=1:p
    x(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end
W = zeros(n*Tt,n*(n-1)/2);
count = 0;
for i=2:n
    W(i:n:end,count+1:count+i-1) = -shortYt(:,1:i-1);
    count = count + i-1;
end
X = SURform2([ones(Tt,1) x],n); 
tmpidx = zeros(k_beta,2);
count = 1;
for i=2:n
    for j=1:i-1
        tmpidx(count,:) = [i j];
        count = count+1;
    end
end
idx_B = sub2ind([n n],tmpidx(:,1),tmpidx(:,2)); 
W = sparse(W);
B = eye(n);
Lfree = [2 3 7 4 8 12];
L = eye(n); 
    
tmpyhat0 = zeros(nsims,2*n+1);  %% [point forecasts, prelike, joint den]
tmpyhat1 = zeros(nsims,2*n+1); 
tmpyhat2 = zeros(nsims,2*n+1); 
tmpyhat4 = zeros(nsims,2*n+1); 

%% initialize
Sigh = .1*ones(n,1);
h0 = log(var(shortYt))';
h = repmat(h0',Tt,1);
beta = zeros(k_beta,1);

for isim = 1:nsims + burnin
  
        % sample alp
    iSig = sparse(1:Tt*n,1:Tt*n,reshape(exp(-h)',Tt*n,1)); 
    XiSig = X'*iSig;
    Kalp = sparse(1:k_alp,1:k_alp,1./Valp) + XiSig*X;    
    CKalp = chol(Kalp,'lower');
    alp_hat = CKalp'\(CKalp\(alp0./Valp + XiSig*(yt-W*beta)));
    alp = alp_hat + CKalp'\randn(k_alp,1); 
        
        % sample beta    
    WiSig = W'*iSig;
    Kbeta = sparse(1:k_beta,1:k_beta,1./Vbeta) + WiSig*W;
    CKbeta = chol(Kbeta,'lower');
    beta_hat = CKbeta'\(CKbeta\(beta0./Vbeta + WiSig*(yt-X*alp)));
    beta = beta_hat + CKbeta'\randn(k_beta,1);    
    
        %% sample h  
    U = reshape(yt-X*alp-W*beta,n,Tt)';    
    for i=1:n
        Ystar = log(U(:,i).^2 + .0001 );
        h(:,i) = SVRW(Ystar,h(:,i),Sigh(i),h0(i));
    end 
    
        %% sample h0
    Kh0 = sparse(1:n,1:n,1./Sigh + 1./Vh);
    h0_hat = Kh0\(ah./Vh + h(1,:)'./Sigh);
    h0 = h0_hat + chol(Kh0,'lower')'\randn(n,1);
    
        %% sample Sigh
    e = h - [h0';h(1:Tt-1,:)];
    Sigh = 1./gamrnd(nuh0+Tt/2,1./(Sh0 + sum(e.^2)'/2));
       
    if isim > burnin
        isave = isim - burnin;
        xtp1 = [1 reshape(shortYt(end:-1:end-p+1,:)',1,n*p)];
        htp1 = h(end,:)'; 
        L(Lfree) = beta;    
        A = reshape(alp,p*n+1,n);
        Atilde = L\A';
        if is_last_miss % if the lastest data are missing, do one more forecast horizon
            htp1 = htp1 + sqrt(Sigh).*randn(n,1);              
            EYtp1 = xtp1*Atilde';
            Sigt = (L\sparse(1:n,1:n,exp(htp1)))/L';
            Ytp1 = EYtp1 + (chol(Sigt,'lower')*randn(n,1))'; 
            xtp1 = [1 Ytp1 xtp1(2:end-n)];
        end
        for tt = 1:5
            htp1 = htp1 + sqrt(Sigh).*randn(n,1);            
            EYtp1 = xtp1*Atilde';
            Sigt = (L\sparse(1:n,1:n,exp(htp1)))/L';
            CSig = chol(Sigt,'lower');
            dSig = diag(Sigt)';
            if tt == 1
                tmpu = CSig\(data_tpk(1,var_small)-EYtp1)';
                lden_joint = -n/2*log(2*pi) - sum(diag(log(CSig)))...
                    -.5*(tmpu'*tmpu);
                lden = -.5*log(2*pi*dSig) -.5*(data_tpk(1,var_small)-EYtp1).^2./dSig; 
                tmpyhat0(isave,:) = [EYtp1 lden lden_joint]; 
            elseif tt == 2 && t<=T-tt
                tmpu = CSig\(data_tpk(2,var_small)-EYtp1)';
                lden_joint = -n/2*log(2*pi) -sum(diag(log(CSig)))...
                    -.5*(tmpu'*tmpu);
                lden = -.5*log(2*pi*dSig) -.5*(data_tpk(2,var_small)-EYtp1).^2./dSig;
                tmpyhat1(isave,:) = [EYtp1 lden lden_joint]; 
            elseif tt == 3 && t<=T-tt
                tmpu = CSig\(data_tpk(3,var_small)-EYtp1)';
                lden_joint = -n/2*log(2*pi) -sum(diag(log(CSig)))...
                    -.5*(tmpu'*tmpu);
                lden = -.5*log(2*pi*dSig) -.5*(data_tpk(3,var_small)-EYtp1).^2./dSig;
                tmpyhat2(isave,:) = [EYtp1 lden lden_joint]; 
            elseif tt == 5 && t<=T-tt
                tmpu = CSig\(data_tpk(5,var_small)-EYtp1)';
                lden_joint = -n/2*log(2*pi) -sum(diag(log(CSig)))...
                    -.5*(tmpu'*tmpu);
                lden = -.5*log(2*pi*dSig) -.5*(data_tpk(5,var_small)-EYtp1).^2./dSig;
                tmpyhat4(isave,:) = [EYtp1 lden lden_joint];
            end         
            Ytp1 = EYtp1 + (chol(Sigt,'lower')*randn(n,1))';
            xtp1 = [1 Ytp1 xtp1(2:end-n)];
        end        
        
    end
end