% This script produces forecasts from a 4-variable BVAR
%
% See:
% Chan, J.C.C. (2020). Large Bayesian Vector Autoregressions. In: P. Fuleky (Eds),
% Macroeconomic Forecasting in the Era of Big Data, 95-125, Springer, Cham

var_small = [1 7 8 12]'; % real GDP growth, IP, unemployment, PCE rate
n = length(var_small);
k = n*p+1;
Y0 = data_t(1:4,var_small);
shortYt = data_t(5:end,var_small);
Yt = reshape(shortYt',Tt*n,1);

tmpyhat0 = zeros(nsims,2*n+1);  %% [point forecasts, prelike, joint den]
tmpyhat1 = zeros(nsims,2*n+1); 

%% construct the Minnesota  prior
[beta_Minn,V_Minn,Sig_hat] = prior_Minn(p,c1,c2,c3,Y0,shortYt);

%% compute a few things
tmpY = [Y0(end-p+1:end,:); shortYt];
Z = zeros(Tt,n*p);
for i=1:p
    Z(:,(i-1)*n+1:i*n) = tmpY(p-i+1:end-i,:);
end
Z = [ones(Tt,1) Z];
X = SURform2(Z,n);
XiSig = X'*kron(speye(Tt),sparse(1:n,1:n,1./Sig_hat));
Kbeta = sparse(1:n*k,1:n*k,1./V_Minn) + XiSig*X;
C_Kbeta = chol(Kbeta,'lower');
beta_hat = C_Kbeta'\(C_Kbeta\(beta_Minn./V_Minn + XiSig * Yt));
CSig = sparse(1:n,1:n,sqrt(Sig_hat));

for isim = 1:nsims + burnin  
        % sample beta    
    beta = beta_hat + C_Kbeta'\randn(k*n,1); 
    
    if isim > burnin
        isave = isim - burnin;
        A = reshape(beta,k,n);
        xtp1 = [1 reshape(shortYt(end:-1:end-p+1,:)',1,n*p)];
        if is_last_miss % if the lastest data are missing, do one more forecast horizon
            Ytp1 = xtp1*A + (CSig*randn(n,1))';
            xtp1 = [1 Ytp1 xtp1(2:end-n)];
        end            
        for tt=1:2
            EYtp1 = xtp1*A;
            dSig = Sig_hat';
            if tt == 1
                tmpu = CSig\(data_tpk(1,var_small)-EYtp1)';
                lden_joint = -n/2*log(2*pi) -sum(diag(log(CSig))) -.5*(tmpu'*tmpu);
                lden = -.5*log(2*pi*dSig) - .5*(data_tpk(1,var_small)-EYtp1).^2./dSig;                
                tmpyhat0(isave,:) = [EYtp1 lden lden_joint]; 
            elseif tt == 2 && t<=T-tt
                tmpu = CSig\(data_tpk(2,var_small)-EYtp1)';
                lden_joint = -n/2*log(2*pi) -sum(diag(log(CSig))) -.5*(tmpu'*tmpu);
                lden = -.5*log(2*pi*dSig) - .5*(data_tpk(2,var_small)-EYtp1).^2./dSig;                
                tmpyhat1(isave,:) = [EYtp1 lden lden_joint];          
            end
            Ytp1 = EYtp1 + (CSig*randn(n,1))';
            xtp1 = [1 Ytp1 xtp1(2:end-n)];
        end        
        
    end
end

