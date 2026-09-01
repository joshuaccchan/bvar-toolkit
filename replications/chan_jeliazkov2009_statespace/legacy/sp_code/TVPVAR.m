% % =======================================================================
% % TVP-VAR example
% % y_t    = mu_t + Gamma_t y_{t-1} + e_t <=> y_t = X_t beta_t +  e_t,
% % beta_t = beta_{t-1} + u_t,
% % e_t ~ N(0,Omega11),  u_t ~ N(0,diag(Omega22)),
% % 
% % Initial conditions: beta_1 ~ N(0,Vbeta).
% % Priors: Omega11 ~ inverse-Wishart(nu01, S01),
% %     each element of Omega22 ~ inverse-gamma(nu02, S02).
% %
% % Data: U.S. GDP growth, inflation, unemployment and interest rate
% %       from 1948:Q1 to 2009:Q4
% % Note: Inflation rate is not annualized
% %
% % See Chan, J.C.C. and Jeliazkov, I. (2009) Efficient Simulation and
% %     Integrated Likelihood Estimation in State Space Models 
% %     
% % Please report any errors to joshuacc.chan@gmail.com
% % =======================================================================

clear; clc;
load USdata.csv;
Y0 = USdata(1:3,:);
shortY = USdata(4:end,:);
[T n] = size(shortY);
Y = reshape(shortY',T*n,1);
p = 1;        % no. of lags
q = n^2*p+n;  % dim of states
Tq = T*q; Tn = T*n; q2 = q^2; nnp1 = n*(n+1);
nloop = 11000; burnin = 1000;
    %% initialize for storage
store_Omega11 = zeros(nloop-burnin,n,n);
store_Omega22 = zeros(nloop-burnin,q);
store_beta = zeros(nloop-burnin,Tq);
    %% initialize the Markov chain
Omega11 = cov(USdata);
Omega22 = .01*ones(q,1); % store only the diag elements
invOmega11 = Omega11\speye(n);
invOmega22 = 1./Omega22;
beta = zeros(Tq,1);
    %% prior
nu01 = n+3; S01 = eye(n);
nu02 = 3;  S02 = .005*ones(q,1);
Vbeta = 5*ones(q,1);
    %% construct and compute a few thigns
X = zeros(T,n*p);
for i=1:p
    X(:,(i-1)*n+1:i*n) = [Y0(3-i+1:end,:); shortY(1:T-i,:)];
end
bigG = SURform([ones(n*T,1) kron(X,ones(n,1))]);
H = speye(Tq,Tq) - sparse(q+1:Tq,1:(T-1)*q,ones(1,(T-1)*q),Tq,Tq);
newnu1 = nu01 + T;
newnu2 = nu02 + (T-1)/2;

disp('Starting MCMC.... ');
disp(' ' );
start_time = clock;
for loop = 1:nloop    
        %% sample beta
    invS = sparse(1:Tq,1:Tq,[1./Vbeta; repmat(invOmega22,T-1,1)]'); 
    K = H'*invS*H;   
    GinvOmega11 = bigG'*kron(speye(T),invOmega11);
    GinvOmega11G = GinvOmega11*bigG;
    invP = K + GinvOmega11G;
    C = chol(invP);                         % so that C'*C = invP
    betahat = C\(C'\(GinvOmega11*Y));
    beta = betahat + C\randn(Tq,1);         % C is upper triangular
        %% sample Omega11
    e1 = reshape(Y-bigG*beta,n,T);
    newS1 = S01 + e1*e1';
    invOmega11 = wishrnd(newS1\speye(n),newnu1);
    Omega11 = invOmega11\speye(n);
        %% sample Omega22
    e2 = reshape(H*beta,q,T);
    newS2 = S02 + sum(e2(:,2:end).^2,2)/2;
    invOmega22 = gamrnd(newnu2,1./newS2);
    Omega22 = 1./invOmega22;  
    if loop>burnin
        i=loop-burnin;
        store_beta(i,:) = beta';
        store_Omega11(i,:,:) = Omega11;
        store_Omega22(i,:) = Omega22';
    end    
    if ( mod( loop, 2000 ) ==0 )
        disp(  [ num2str( loop ) ' loops... ' ] )
    end    
end
disp( ['MCMC takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
disp(' ' );

    %% compute posterior means
betahat = mean(store_beta)';
Omega11hat = squeeze(mean(store_Omega11)); 
Omega22hat = mean(store_Omega22); 
tempbetahat = reshape(betahat,nnp1,T)';
gamhat1 = tempbetahat(:,1:n+1);
gamhat2 = tempbetahat(:,(n+1)+1:2*(n+1));
gamhat3 = tempbetahat(:,2*(n+1)+1:3*(n+1));
gamhat4 = tempbetahat(:,3*(n+1)+1:end);

    %% plot graphs
figure;
set(0,'DefaultAxesColorOrder',[0 0 0], 'DefaultAxesLineStyleOrder','-.*|-|-.|--|:')
set(0,'DefaultAxesColorOrder',[0 0 1;0 0 0;1 0 0; 0 1 0; 1 0 1])
t = (1948.75:.25:2009.75)';
subplot(2,2,1); 
    plot(t,gamhat1(:,1),'-.',t,gamhat1(:,2),'-',...
        t,gamhat1(:,3),'-.',t,gamhat1(:,4),'--',...
        t,gamhat1(:,5),':','Linewidth',2,'MarkerSize',3);
    title('GDP growth equation'); box off; xlim([1947 2011]);
subplot(2,2,4); plot(t,gamhat2(:,1),'-.',t,gamhat2(:,2),'-',...
    t,gamhat2(:,3),'-.',t,gamhat2(:,4),'--',...
    t,gamhat2(:,5),':','Linewidth',2,'MarkerSize',3);
    title('Inflation rate equation'); box off; xlim([1947 2011]);
subplot(2,2,2); plot(t,gamhat3(:,1),'-.',t,gamhat3(:,2),'-',...
    t,gamhat3(:,3),'-.',t,gamhat3(:,4),'--',...
    t,gamhat3(:,5),':','Linewidth',2,'MarkerSize',3);
    title('Unemployment equation'); box off; xlim([1947 2011]);
subplot(2,2,3); plot(t,gamhat4(:,1),'-.',t,gamhat4(:,2),'-',...
    t,gamhat4(:,3),'-.',t,gamhat4(:,4),'--',...
    t,gamhat4(:,5),':','Linewidth',2,'MarkerSize',3);
    title('Interest rate equation');  box off; xlim([1947 2011]);
h = legend('\mu_{t-1}','GDP_{t-1}','inflation_{t-1}','unemp_{t-1}', ...
    'interest_{t-1}','Location','SouthOutside','Orientation','horizontal');


