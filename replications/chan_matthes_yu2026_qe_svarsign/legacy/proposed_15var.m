% This script produces the empirical application in Chan (2022) 
% using the proposed algorithm: estimating a 15-variable VAR 
% identified with sign restrictions

p = 5;          % if p > 8, need to change Y0 and Y below
nsim = 1000;    % # of posterior draws (that satisfy all the restrictions)
dataset = 2;    % 1: 6-variable; 2: 15-variable
addpath('./data');
addpath('./utility');
nbatch = 1000; % # of posterior draws sampled in a batch
horizon = 36;   % # of steps for impulse responses
    % load data
data = readmatrix('./data/database_2019Q4.csv','Range','B3:P142'); % 15 variables + 5 shocks % shocks: supply, demand, monetary, investment, financial
var_id = 1:15;
idx_ns = [1,2,4,5,10,11,12,13,15]; % index for variables in levels
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
elseif dataset == 3
        % find the optimal kappa values
    [ml_opt,kappa] = get_OptKappa(Y0,Y,Z,p,[.04,.0016],'redu',idx_ns); 
end
sig2 = get_resid_var(Y0,Y);
prior_stru = prior_ACP_stru(n,p,kappa,sig2,idx_ns);
prior_redu = prior_ACP_redu(n,p,kappa,sig2,idx_ns);
% setup the sign restrictions and row inequalities
        % sign restrictions
    supply = [1,-1,NaN,NaN,1,NaN,NaN,NaN,NaN,-1,-1,NaN,1,NaN,1]';
    demand = [1,1,1,NaN,NaN,NaN,NaN,NaN,NaN,1,1,NaN,1,1,NaN]';
    monetary = [1,1,-1,NaN,NaN,NaN,NaN,NaN,NaN,1,1,NaN,1,-1,NaN]';
    invest = [1,1,1,NaN,-1,NaN,NaN,NaN,NaN,1,1,NaN,1,1,-1]';
    finc = [1,1,1,NaN,1,NaN,NaN,NaN,NaN,1,1,NaN,1,1,1]';
    S = [supply,demand,monetary,invest,finc];
    m = size(S,2); % # of shocks
        % row inequalities    
    Rineq = [zeros(1,n); -1,0,0,1,zeros(1,11); zeros(1,n); ...
        1,0,0,-1,zeros(1,11);1,0,0,-1,zeros(1,11)];    
    k = 1;
	
% estimation
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
                if k == 1
                    if (sum(signL(:,j) == S(idx,i)) == nidx) && ...
                            (sum(Rineq(i,:,:)*L(:,j) <= 0) == 1)
                        satTab(i,j) = 1;
                    elseif (sum(signL(:,j) == -S(idx,i)) == nidx) && ...
                            (sum(Rineq(i,:,:)*(-L(:,j)) <= 0) == 1)
                        satTab(i,j) = -1;
                    end
                elseif k>1
                    if (sum(signL(:,j) == S(idx,i)) == nidx) && ...
                            (sum(squeeze(Rineq(i,:,:))'*L(:,j) <= 0) == k)
                        satTab(i,j) = 1;
                    elseif (sum(signL(:,j) == -S(idx,i)) == nidx) && ...
                            (sum(squeeze(Rineq(i,:,:))'*(-L(:,j)) <= 0) == k)
                        satTab(i,j) = -1;
                    end
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

