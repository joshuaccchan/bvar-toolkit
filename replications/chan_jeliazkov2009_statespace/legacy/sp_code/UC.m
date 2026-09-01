% % =======================================================================
% % Unobserved components model with stochastic volatility 
% %
% % y_t    = tau_t + e_t,      e_t ~ N(0,exp(h_t)),
% % tau_t  = tau_{t-1} + u_t,  u_t ~ N(0,omega2tau), 
% % h_t    = h_{t-1} + v_t,    v_t ~ N(0,omega2h).
% % 
% % Initial conditions: tau_1 ~ N(0,Vtau), h_1 ~ N(0,Vh).
% % Priors: omega2tau ~ inverse-gamma(nutau, Stau),
% %         omega2h   ~ inverse-gamma(nuh, Sh).
% % Data: U.S. CPI inflation from 1947:Q2 to 2011:Q3
% %
% % See Chan, J.C.C. (2012) Moving average stochastic volatility models
% %     with application to inflation forecast
% %
% % Please report any errors to joshuacc.chan@gmail.com
% % =======================================================================

clear; clc;
nloop = 11000;
burnin = 1000;
load 'USCPI.csv';
y = USCPI; 
T = length(y);
    %% prior
Vtau = 9; Vh = 9;
nutau = 5; Stau = .25^2*(nutau-1);
nuh = 5; Sh = .2^2*(nuh-1);
    %% initialize the Markov chain
omega2tau = .25^2;
omega2h = .2^2;
h = log(var(y)*.8)*ones(T,1);
    %% initialize for storeage
store_omega2tau = zeros(nloop - burnin,1); 
store_omega2h = zeros(nloop - burnin,1); 
store_tau = zeros(nloop - burnin,T); 
store_h = zeros(nloop - burnin,T);
    %% compute a few things outside the loop
H = speye(T) - sparse(2:T,1:(T-1),ones(1,T-1),T,T);    
newnutau = (T-1)/2 + nutau;
newnuh = (T-1)/2 + nuh;

disp('Starting MCMC.... ');
disp(' ' );
start_time = clock;
for loop = 1:nloop
        %% sample tau    
    invOmegatau = sparse(1:T,1:T,[1/Vtau 1/omega2tau*ones(1,T-1)]);
    invSigy = sparse(1:T,1:T,exp(-h));
    Ktau = H'*invOmegatau*H + invSigy;
    Ctau = chol(Ktau,'lower');          % so that Ctau*Ctau' = Ktau
    tauhat = Ktau\(invSigy*y);
    tau = tauhat + Ctau'\randn(T,1);    % note the transpose 
        %% sample h
    ystar = log((y-tau).^2 + .0001 );
    h = SVRW(ystar,h,omega2h,Vh);   
        %% sample omega2tau
    newStau = Stau + sum((tau(2:end)-tau(1:end-1)).^2)/2;
    omega2tau = 1/gamrnd(newnutau, 1./newStau);        
        %% sample omega2h
    newSh = Sh + sum((h(2:end)-h(1:end-1)).^2)/2;
    omega2h = 1/gamrnd(newnuh, 1./newSh);    
    if loop>burnin
        i = loop-burnin;
        store_tau(i,:) = tau';
        store_h(i,:) = h'; 
        store_omega2tau(i,:) = omega2tau;
        store_omega2h(i,:) = omega2h;    
    end    
    if ( mod( loop, 2000 ) ==0 )
        disp(  [ num2str( loop ) ' loops... ' ] )
    end    
end
disp( ['MCMC takes '  num2str( etime( clock, start_time) ) ' seconds' ] );
disp(' ' );

tauhat = mean(store_tau)';
hhat = mean(store_h)'; 

tid = (1947.25:.25:2011.5)';   
figure;    
subplot(1,2,1); plot(tid, tauhat, 'LineWidth',2,'Color','black');        
xlim([1947 2012]);  title('\tau_t');
subplot(1,2,2); plot(tid, hhat, 'LineWidth',2,'Color','black');        
xlim([1947 2012]);   title('h_t');
set(gcf,'Position',[100 100 800 300]);