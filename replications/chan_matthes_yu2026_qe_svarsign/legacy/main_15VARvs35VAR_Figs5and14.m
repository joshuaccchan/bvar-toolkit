% This is the main script for replicating the application in section 4 (35-VAR vs. 15-VAR): Figures 5 and 14 (in Appendix B)
% NOTICE: To run the code successfully, one needs to first execute main_15VAR_Fig1.m to generate the output file results_15var.mat,
% which is required as input for the subsequent 15-VAR analysis.

clear; clc;
rng(123);
p = 5;          % if p > 8, need to change Y0 and Y below
dataset = 1;    % 1: 35-variable
nsim = 1000;    % # of posterior draws (that satisfy all the restrictions)

addpath('./data');
addpath('./utility');
nbatch = 1000; % # of posterior draws sampled in a batch
horizon = 36;   % # of steps for impulse responses
% load data
if dataset == 1
    data = readmatrix('./data/data35_1975.csv','Range','B35:AJ182');
    % 35 variables + 8 shocks
    % shocks: Demand, Investment, Financial, Monetary, Govt spending
    %   Technology, Labor supply, Wage bargaining
    var_id = 1:35;
    idx_ns = 1:35; % index for variables in levels
end
data = data(9:end,:); %1985Q1:2019Q4, consistent with the 15-var
data(:,32) = data(:,32) - data(:,25); %create spread, which is the BAA-FFR
data(:,29) = data(:,32) - data(:,29); %create spread2, which is the BAA-10yearTbillrate
data_SP500 = readmatrix('./data/database_2019Q4.csv','Range','F3:F142');
data(:,34) = data_SP500; %replace SP500 with real SP500 to match the 15 var case
% and we use PNFIx (nonresidential investment) to surrogate the investment in the 15-var
Y0 = data(1:8,:);  % save the first 8 obs as the initial conditions
Y = data(9:end,:);
[T,n] = size(Y);
tmpY = [Y0(end-p+1:end,:); Y];
Z = zeros(T,n*p);
for ii=1:p
    Z(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
end
Z = [ones(T,1) Z];
if dataset == 1
    % find the optimal kappa values
    [ml_opt,kappa] = get_OptKappa(Y0,Y,Z,p,[.04,.0016],'redu',idx_ns);
    % [ml_opt,kappa] = get_OptKappa_ver2(Y0,Y,Z,p,[.04,.0016,1],'redu',idx_ns);
end
sig2 = get_resid_var(Y0,Y);
prior_stru = prior_ACP_stru(n,p,kappa,sig2,idx_ns);
prior_redu = prior_ACP_redu(n,p,kappa,sig2,idx_ns);

% setup the sign restrictions and row inequalities
if dataset == 1
    % sign restrictions
    demand = [1,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,1,1,1,1,NaN,NaN,NaN,NaN,-1,...
        1,1,NaN,NaN,NaN,1,1,NaN,NaN,NaN,1,NaN,NaN,NaN,NaN,NaN]';
    investment = [1,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,1,1,1,1,NaN,NaN,NaN,NaN,-1,...
        1,1,NaN,NaN,NaN,1,1,NaN,NaN,NaN,1,NaN,NaN,NaN,-1,NaN]';
    financial = [1,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,1,1,1,1,NaN,NaN,NaN,NaN,-1,...
        1,1,NaN,NaN,NaN,1,1,NaN,NaN,NaN,1,NaN,NaN,NaN,1,NaN]';
    monetary = [-1,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,-1,-1,-1,-1,-1,NaN,NaN,NaN...
        -1,1,-1,-1,NaN,NaN,NaN,1,1,1,1,1,1,1,NaN,NaN,-1,NaN]';
    government = [1,NaN,NaN,NaN,NaN,NaN,1,-1,1,1,1,1,1,1,NaN,NaN,NaN,NaN,-1,...
        NaN,NaN,NaN,NaN,NaN,1,1,NaN,NaN,NaN,1,NaN,NaN,NaN,NaN,NaN]';
    technology = [1,1,NaN,1,NaN,NaN,NaN,NaN,NaN,-1,-1,-1,-1,-1,1,1,1,NaN,-1,...
        NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN]';
    labor = [1,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,-1,-1,-1,-1,-1,-1,NaN,NaN,NaN,1,...
        NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN]';
    wage = [1,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,-1,-1,-1,-1,-1,-1,NaN,NaN,NaN,-1,...
        NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN]';
    S = [demand,investment,financial,monetary,government,technology,labor,wage];
    m = size(S,2); % # of shocks
    % row inequalities
    k = 2; % # of different types of ranking restrictions
    Rineq = zeros(m,n,2);
    % ranking of investment and output
    Rineq(1,1,1) = -1; Rineq(1,4,1) = 1; % demand shock: log_investment - log_output < 0
    Rineq(2,1,1) = 1; Rineq(2,4,1) = -1; % investment shock: -log_investment + log_output < 0
    Rineq(3,1,1) = 1; Rineq(3,4,1) = -1; % financial shock: -log_investment + log_output < 0
    % ranking of government spending and output
    Rineq(1,1,2) = -1; Rineq(1,7,2) = 1; % demand shock: log_government - log_output < 0
    Rineq(2,1,2) = -1; Rineq(2,7,2) = 1; % investment shock: log_government - log_output < 0
    Rineq(3,1,2) = -1; Rineq(3,7,2) = 1; % financial shock: log_government - log_output < 0
    Rineq(5,1,2) = 1; Rineq(5,7,2) = -1; % G spending shock: -log_government + log_output < 0
end

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

        [Q,R] = QR(randn(n,n));
        L0 = chol(Sigtilde,'lower');
        L = L0*Q;
        satTab = zeros(m,n); % (i,j) = 1 if the j-th column of L satisfies all restrictions for the i-th shock
        % (i,j) = -1 if the negative of j-th column of L satisfies all restrictions
        for i=1:m  % check sign restrictions & row inequilities
            idx = find(S(:,i)==-1 | S(:,i)==1);
            nidx = length(idx);
            signL = sign(L(idx,:));
            for j=1:n
                if (sum(signL(:,j) == S(idx,i)) == nidx) && ...
                        (sum(squeeze(Rineq(i,:,:))'*L(:,j) <= 0) == k)
                    satTab(i,j) = 1;
                elseif (sum(signL(:,j) == -S(idx,i)) == nidx) && ...
                        (sum(squeeze(Rineq(i,:,:))'*(-L(:,j)) <= 0) ==k)
                    satTab(i,j) = -1;
                end

            end
        end
        if nnz(sum(abs(satTab),2)) == m  % admissible set is non-empty
            reorder = zeros(1,n);
            for i=1:m
                idx = find(satTab(i,:));
                draw = idx(unidrnd(length(idx)));
                reorder(i) = draw;
                if satTab(i,draw) == -1
                    L(:,draw) = -L(:,draw);
                end
            end
            reorder(m+1:end) = setdiff(1:n,reorder);
            L = L(:,reorder);
            % randomly permute and switch the signs of the last n-m columns
            tmpL2 = L(:,m+1:end);
            tmpL2 = tmpL2(:,randperm(n-m))*diag(2*(rand(n-m,1)>.5) - 1);
            L(:,m+1:end) = tmpL2;
            count_sat = count_sat + 1;
            Btilde = reshape(store_Btilde(isim,:),n*p+1,n);
            response = IRredu(Btilde(2:end,:),L,horizon,m);
            store_response(count_sat,:,:,:) = response;
        end
    end

    if (mod(count_total, 1000) == 0)
        disp(['Out of ' num2str(count_total) ' posterior draws, ' ...
            num2str(count_sat) ' satisfy all the restrictions.']);
    end
end
disp([num2str(nsim) ' posterior draws that satisfy all the restrictions are obtained.']);
disp(['The simulation took ' num2str(etime(clock,start_time)/60) ' minutes.']);

response_median_35 = squeeze(median(store_response));
response_CI_35 = squeeze(quantile(store_response,[.16,.84]));

save('results_35var_compare15var.mat','response_median_35','response_CI_35')

load('results_15var.mat')
response_median_15 = squeeze(median(store_response));
response_CI_15 = squeeze(quantile(store_response,[.16,.84]));

IRF_varid = [1,10,26,4,34,32];
titletext = ["GDP" "GDP Deflator" "3-month Tbill" "Investment" "S&P 500" "Spread"];
shocktext = ["Demand" "Investment"  "Financial" "Monetary" ...
    "Labor Supply" ];
scaleno = [0.01, 0.01, 1, 0.01, 1, 1];
for shock_l = 3
    figure;
    for jj = 1:6
        subplot(3,2,jj);
        hold on
        plotCI((0:horizon-1)',squeeze(response_CI_35(1,IRF_varid(jj),shock_l,:))*scaleno(jj),...
            squeeze(response_CI_35(2,IRF_varid(jj),shock_l,:))*scaleno(jj));
        plotCI_othercolor((0:horizon-1)',squeeze(response_CI_15(1,jj,5,:)),...
            squeeze(response_CI_15(2,jj,5,:)),'r');
        h2 = plot(0:horizon-1,squeeze(response_median_35(IRF_varid(jj),shock_l,:))*scaleno(jj),'--k','LineWidth',1, 'DisplayName','35-VAR');
        h1 = plot(0:horizon-1,squeeze(response_median_15(jj,5,1:end)),'--r','LineWidth',1, 'DisplayName', '15-VAR         ');

        hold off
        line(xlim, [0,0], 'Color', 'k', 'LineWidth', .5); % Draw line for x-axis
        xlim([-.1, horizon-1]);
        if jj == 1 ylim([-0.001, 0.008]); end
        if jj == 2 ylim([-0.001, 0.003]); end
        if jj == 3 ylim([-0.2, 0.4]); end
        if jj == 4 ylim([-0.01, 0.03]); end
        if jj == 5 ylim([-0.005, 0.015]); end
        if jj == 6 ylim([-0.4, 0.2]); end
        grid on; title(titletext(jj));
        box off;
    end
    if jj== 6
        legend([h1, h2], 'Location', 'southoutside', 'Orientation', 'horizontal');
        set(legend, 'Position', [0.4, 0.01, 0.2, 0.05]);
    end
    sgtitle(shocktext(shock_l)) ;
    set(gcf,'Position',[100 300 600 400]);
end
print(gcf, 'Fig5', '-depsc', '-painters');

for shock_l = 3
    figure;
    for jj = 1:6
        hold on
        len1 = squeeze(response_CI_15(2,jj,5,:))-squeeze(response_CI_15(1,jj,5,:));
        len2 = (squeeze(response_CI_35(2,IRF_varid(jj),shock_l,:))-squeeze(response_CI_35(1,IRF_varid(jj),shock_l,:)))*scaleno(jj);
        ratio = len2./len1;
        plot(0:horizon-1,ratio, 'LineWidth', 1);
        hold off

        xlim([-.1, horizon-1]);
        grid on;
        box off;
    end

end
h = yline(1, '--', 'LineWidth', 1, 'Color', 'k');
h.HandleVisibility = 'off';
legend(["GDP" "GDP Deflator" "3-month Tbill" "Investment" "S&P 500" "Spread"], ...
    'Location', [0.15, 0.78, 0.21, 0.00]);
sgtitle(shocktext(shock_l)) ;
set(gcf,'Position',[100 300 600 400]);
print(gcf, 'Fig14', '-depsc', '-painters');

