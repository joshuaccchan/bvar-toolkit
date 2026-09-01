% This script reproduces the empirical application in Chan (2022):
% estimating a 6- or 15-variable VAR identified with sign restrictions
%
% This code is free to use for academic purposes only, provided that the 
% paper is cited as:
%
% Chan, J.C.C. (2022). Asymmetric Conjugate Priors for Large Bayesian VARs,
% Quantitative Economics, 13(3): 1145-1169
%
% This code comes without technical support of any kind. It is expected to
% reproduce the results reported in the paper. Under no circumstances will
% the authors be held responsible for any use (or misuse) of this code in
% any way.

clear; clc;
p = 5;          % if p > 8, need to change Y0 and Y below
dataset = 1;    % 1: 6-variable; 2: 15-variable
if dataset == 1
    nsim = 5000;   % # of posterior draws (that satisfy all the restrictions)
elseif dataset == 2
    nsim = 1000;    % might take a few days to get 1000 draws
end

addpath('./utility');
nbatch = 50000; % # of posterior draws sampled in a batch
horizon = 36;   % # of steps for impulse responses
    % load data
data = xlsread('database_2019Q4.xlsx');
if dataset == 1
        % 6 variables + 5 shocks
        % shocks: supply, demand, monetary, investment, financial
    var_id = 1:6; %  GDP, deflator, interest rate, investiment, S&P, spread 1    
    idx_ns = [1,2,4,5]; % index for variables in levels
elseif dataset == 2
        % 15 variables + 5 shocks
        % shocks: supply, demand, monetary, investment, financial
    var_id = 1:15;
    idx_ns = [1,2,4,5,10,11,12,13,15]; % index for variables in levels
end
Y0 = data(1:8,var_id);  % save the first 8 obs as the initial conditions
Y = data(9:end,var_id);
[T,n] = size(Y);
tmpY = [Y0(end-p+1:end,:); Y];
Z = zeros(T,n*p); 
for ii=1:p
    Z(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
end
Z = [ones(T,1) Z];
if dataset == 1
    kappa = [1, 1, 1, 100]; % [kappa1, kappa2, kappa3, kappa4]
elseif dataset == 2
        % find the optimal kappa values
    [ml_opt,kappa] = get_OptKappa(Y0,Y,Z,p,[.04,.0016],'redu',idx_ns); 
end
sig2 = get_resid_var(Y0,Y);
prior_stru = prior_ACP_stru(n,p,kappa,sig2,idx_ns);
prior_redu = prior_ACP_redu(n,p,kappa,sig2,idx_ns);

    % setup the sign restrictions and row inequalities
if dataset == 1
        % sign restrictions
    supply = [1,-1,NaN,NaN,1,NaN]';
    demand = [1,1,1,NaN,NaN,NaN]';
    monetary = [1,1,-1,NaN,NaN,NaN]';
    invest = [1,1,1,NaN,-1,NaN]';
    finc = [1,1,1,NaN,1,NaN]';
    S = [supply,demand,monetary,invest,finc];
    m = size(S,2); % # of shocks
        % row inequalities    
    Rineq = [-1,0,0,1,0,0;1,0,0,-1,0,0;1,0,0,-1,0,0];
    Ridx = [2,4,5];    % column indices Rineq applies to 
elseif dataset == 2
        % sign restrictions
    supply = [1,-1,NaN,NaN,1,NaN,NaN,NaN,NaN,-1,-1,NaN,1,NaN,1]';
    demand = [1,1,1,NaN,NaN,NaN,NaN,NaN,NaN,1,1,NaN,1,1,NaN]';
    monetary = [1,1,-1,NaN,NaN,NaN,NaN,NaN,NaN,1,1,NaN,1,-1,NaN]';
    invest = [1,1,1,NaN,-1,NaN,NaN,NaN,NaN,1,1,NaN,1,1,-1]';
    finc = [1,1,1,NaN,1,NaN,NaN,NaN,NaN,1,1,NaN,1,1,1]';
    S = [supply,demand,monetary,invest,finc];
    m = size(S,2); % # of shocks
        % row inequalities    
    Rineq = [-1,0,0,1,zeros(1,11);1,0,0,-1,zeros(1,11);1,0,0,-1,zeros(1,11)]; 
    Ridx = [2,4,5];    % column indices Rineq applies to 
end
nR = length(Ridx); % # of row inequalities 

start_time = clock;
count_sat = 0; % counter for # draws that satisfy all the conditions
count_total = 0; % counter for total # draws
disp(['Computing impulse responses from a ' num2str(n) '-variable VAR']);
disp('    to an one-standard-deviation financial shock...');
store_response = zeros(nsim,n,m,horizon);
while count_sat < nsim     
        %  sample nbatch draws from the posterior
    [store_alp,store_beta,store_Sig] = sample_ThetaSig(Y0,Y,p,prior_redu,nbatch);
    count_total = count_total + nbatch;
    
        % obtain the reduced-form parameters
    [store_Btilde,store_Sigtilde] = getReducedForm(store_alp,store_beta,store_Sig);    
    
    for isim = 1:nbatch  % go through the nbatch posterior draws to find those 
        Sigtilde = squeeze(store_Sigtilde(isim,:,:));        
        msat = 0;  % counter for the # of shocks that satisfies the sign restrictions
        nRsat = 0; % counter for the # of satisfied row inequalities
        
        [Q,R] = QR(randn(n,n));
        L0 = chol(Sigtilde,'lower');
        L = L0*Q;
        for i=1:m  % check sign restrictions
            idx = find(S(:,i)==-1 | S(:,i)==1);
            nidx = length(idx);
            signL = sign(L(idx,:));
                % check if the i-th column satisfies the sign restrictions
            if (sum(signL(:,i) == S(idx,i)) == nidx) 
                msat = msat + 1;
                % or if the negative of the i-th column satisfies the sign restrictions
            elseif (sum(signL(:,i) == -S(idx,i)) == nidx)
                L(:,i) = -L(:,i); % change the sign of the i-th column
                msat = msat + 1;
            else
                break
            end
        end
        for j=1:nR % check row inequalities
            if Rineq(j,:)*L(:,Ridx(j)) < 0
                nRsat = nRsat + 1;
            else
                nRsat = 0;
                break
            end
        end
        if msat == m && nRsat == nR            
            count_sat = count_sat + 1;
            Btilde = reshape(store_Btilde(isim,:),n*p+1,n);
            response = IRredu(Btilde(2:end,:),L,horizon,m);
            store_response(count_sat,:,:,:) = response;
        end
    end
    
    if (mod(count_total, 10^6) == 0)
        disp(['Out of ' num2str(count_total/10^6) ' million posterior draws, ' ...
            num2str(count_sat) ' satisfy all the restrictions.']);
    end   
end
disp([num2str(nsim) ' posterior draws that satisfy all the restrictions are obtained.']);
disp(['The simulation took ' num2str(etime(clock,start_time)/60) ' minutes.']);

response_median = squeeze(median(store_response));
response_CI = squeeze(quantile(store_response,[.16,.84]));
ylim_u = [.008, .003, .4, .03, .015, .2];
ylim_l = [-.001, -.001, -.2, -.01, -.005, -.4];
if dataset == 1
    titletext = ["GDP" "GDP Deflator" "3-month Tbill" "Investment" "S&P 500" "Spread"];
    figure;
    for jj = 1:6
        subplot(3,2,jj);
        hold on
        plotCI((1:horizon-1)',squeeze(response_CI(1,jj,5,2:end)),squeeze(response_CI(2,jj,5,2:end)));
        plot(1:horizon-1,squeeze(response_median(jj,5,2:end)),'--k','LineWidth',1);
        hold off
        line(xlim, [0,0], 'Color', 'k', 'LineWidth', .5); % Draw line for x-axis
        grid on; title(titletext(jj)); xlim([.95, horizon-1]); 
        ylim([ylim_l(jj) ylim_u(jj)]); box off;
    end
    set(gcf,'Position',[100 300 600 400]);
elseif dataset == 2
    titletext = ["GDP" "GDP Deflator" "3-month Tbill" "Investment" "S&P 500" "Spread"];    
    figure;
    for jj = 1:6
        subplot(3,2,jj);
        hold on
        plotCI((1:horizon-1)',squeeze(response_CI(1,jj,5,2:end)),squeeze(response_CI(2,jj,5,2:end)));
        plot(1:horizon-1,squeeze(response_median(jj,5,2:end)),'--k','LineWidth',1);
        hold off
        line(xlim, [0,0], 'Color', 'k', 'LineWidth', .5); % Draw line for x-axis
        grid on; title(titletext(jj)); xlim([.95, horizon-1]); 
        ylim([ylim_l(jj) ylim_u(jj)]); box off;
    end    
    set(gcf,'Position',[100 300 600 400]);
    
    titletext = ["GDP Deflator" "CPI" "PCE"];    
    varidx = [2,10,11];
    figure;
    for jj = 1:3
        subplot(1,3,jj);
        hold on
        plotCI((1:horizon-1)',squeeze(response_CI(1,varidx(jj),5,2:end)),...
            squeeze(response_CI(2,varidx(jj),5,2:end)));
        plot(1:horizon-1,squeeze(response_median(varidx(jj),5,2:end)),'--k','LineWidth',1);
        hold off
        line(xlim, [0,0], 'Color', 'k', 'LineWidth', .5); % Draw line for x-axis
        ylim([-.001 .003]); box off;
        grid on; title(titletext(jj)); xlim([.95, horizon-1]); box off;
    end   
    set(gcf,'Position',[100 100 600 200])
end