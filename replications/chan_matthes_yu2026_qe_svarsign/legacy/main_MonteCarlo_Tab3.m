% This is the main script for replicating the Monte Carlo experiments in section 3: Table 3
% NOTICE: Full execution is time-intensive. Setting time_set = 300 yields approximate results,
% which are subject to approximation and Monte Carlo errors.

clear; clc;
rng(123);
addpath('./utility');
addpath('./Restr');
time_set = 300; % time limit for each case, default is 10000 seconds
p = 5; nbatch = 100; % # of posterior draws sampled in a batch

%%% ---------------------------------- Panel A ----------------------------------%%%
Tab3_PanelA = zeros(4,6);
for i_case = 1:6
    eval(['load', ' SignRestr', num2str(i_case), '.mat']); eval(['load', ' DataT200_', num2str(i_case), '.mat']); % load data and restriction matrix

    [T, n] = size(data); var_id = 1:n;
    Y0 = data(1:8,var_id); Y = data(9:end,var_id);
    [T,n] = size(Y); tmpY = [Y0(end-p+1:end,:); Y];
    Z = zeros(T,n*p);
    for ii=1:p
        Z(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
    end
    Z = [ones(T,1) Z]; kappa = [1, 1, 1, 100]; sig2 = get_resid_var(Y0,Y);
    prior_stru = prior_ACP_stru(n,p,kappa,sig2,[]); prior_redu = prior_ACP_redu(n,p,kappa,sig2,[]);

    S = Restr; m = size(S,2); % setup the sign restrictions
    Rineq = zeros(m,n); Ridx = [1:m]; nR = length(Ridx); % setup the row inequalities
    %%% ------------------RWZ method-------------%%%
    time_cost = tic; count_sat = 0; count_total = 0;
    while toc(time_cost) < time_set
        [store_alp,store_beta,store_Sig] = sample_ThetaSig(Y0,Y,p,prior_redu,nbatch);
        count_total = count_total + nbatch;
        [store_Btilde,store_Sigtilde] = getReducedForm(store_alp,store_beta,store_Sig);
        for isim = 1:nbatch
            Sigtilde = squeeze(store_Sigtilde(isim,:,:)); msat = 0; nRsat = 0;
            [Q,R] = QR(randn(n,n)); L0 = chol(Sigtilde,'lower'); L = L0*Q;
            for i=1:m  % check sign restrictions
                idx = find(S(:,i)==-1 | S(:,i)==1); nidx = length(idx); signL = sign(L(idx,:));
                if (sum(signL(:,i) == S(idx,i)) == nidx)
                    msat = msat + 1;
                elseif (sum(signL(:,i) == -S(idx,i)) == nidx)
                    L(:,i) = -L(:,i);
                    msat = msat + 1;
                else
                    break
                end
            end
            for j=1:nR % check row inequilities
                if Rineq(j,:)*L(:,Ridx(j)) <= 0
                    nRsat = nRsat + 1;
                else
                    nRsat = 0;
                    break
                end
            end
            if msat == m && nRsat == nR
                count_sat = count_sat + 1;
            end
        end
        if (mod(count_total, 10000) == 0)
            disp(['Out of ' num2str(count_total) ' posterior draws, ' num2str(count_sat) ' satisfy all the restrictions.']);
        end
    end
    Results_RWZ = [count_total/10^6, count_sat]*10000/time_set
    %%% ------------------proposed method-------------%%%
    time_cost = tic; count_sat = 0; count_total = 0; k = 1;
    while toc(time_cost) < time_set
        [store_alp,store_beta,store_Sig] = sample_ThetaSig(Y0,Y,p,prior_redu,nbatch);
        count_total = count_total + nbatch;
        [store_Btilde,store_Sigtilde] = getReducedForm(store_alp,store_beta,store_Sig);
        for isim = 1:nbatch
            Sigtilde = squeeze(store_Sigtilde(isim,:,:));
            [Q,R] = QR(randn(n,n)); L0 = chol(Sigtilde,'lower'); L = L0*Q; satTab = zeros(m,n);
            for i=1:m  % check sign restrictions & row inequilities
                idx = find(S(:,i)==-1 | S(:,i)==1); nidx = length(idx); signL = sign(L(idx,:));
                for j=1:n
                    if (sum(signL(:,j) == S(idx,i)) == nidx) && ...
                            (sum((Rineq(i,:))*L(:,j) <= 0) == k) 
                        satTab(i,j) = 1;
                    elseif (sum(signL(:,j) == -S(idx,i)) == nidx) && ...
                            (sum((Rineq(i,:))*(-L(:,j)) <= 0) ==k) 
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
                reorder(m+1:end) = setdiff(1:n,reorder); L = L(:,reorder); % randomly permute and switch the signs of the last n-m columns
                tmpL2 = L(:,m+1:end); tmpL2 = tmpL2(:,randperm(n-m))*diag(2*(rand(n-m,1)>.5) - 1); L(:,m+1:end) = tmpL2;
                count_sat = count_sat + 1;
            end
        end
    end
    Results_proposed = [count_total/10^6, count_sat]*10000/time_set

    Tab3_PanelA(:,i_case) = [Results_RWZ, Results_proposed]
end
Tab3_PanelA = [Tab3_PanelA(:,1:3); Tab3_PanelA(:,4:6)]

%%% ---------------------------------- Panel B ----------------------------------%%%
Tab3_PanelB = zeros(4,6);
for i_case = 1:6
    eval(['load', ' SignRestr', num2str(i_case), '.mat']); eval(['load', ' DataT200_', num2str(i_case), '.mat']); % load data and sign restriction matrix
    eval(['load', ' RankRestr', num2str(i_case), '.mat']); % load ranking restriction matrix

    [T, n] = size(data); var_id = 1:n;
    Y0 = data(1:8,var_id); Y = data(9:end,var_id);
    [T,n] = size(Y); tmpY = [Y0(end-p+1:end,:); Y];
    Z = zeros(T,n*p);
    for ii=1:p
        Z(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
    end
    Z = [ones(T,1) Z]; kappa = [1, 1, 1, 100]; sig2 = get_resid_var(Y0,Y);
    prior_stru = prior_ACP_stru(n,p,kappa,sig2,[]); prior_redu = prior_ACP_redu(n,p,kappa,sig2,[]);

    S = Restr; m = size(S,2); % setup the sign restrictions
    %%% ------------------RWZ method-------------%%%
	Rineq = Rineq_RWZ; nR = length(Ridx); % setup the row inequalities
    time_cost = tic; count_sat = 0; count_total = 0;
    while toc(time_cost) < time_set
        [store_alp,store_beta,store_Sig] = sample_ThetaSig(Y0,Y,p,prior_redu,nbatch);
        count_total = count_total + nbatch;
        [store_Btilde,store_Sigtilde] = getReducedForm(store_alp,store_beta,store_Sig);
        for isim = 1:nbatch
            Sigtilde = squeeze(store_Sigtilde(isim,:,:)); msat = 0; nRsat = 0;
            [Q,R] = QR(randn(n,n)); L0 = chol(Sigtilde,'lower'); L = L0*Q;
            for i=1:m  % check sign restrictions
                idx = find(S(:,i)==-1 | S(:,i)==1); nidx = length(idx); signL = sign(L(idx,:));
                if (sum(signL(:,i) == S(idx,i)) == nidx)
                    msat = msat + 1;
                elseif (sum(signL(:,i) == -S(idx,i)) == nidx)
                    L(:,i) = -L(:,i);
                    msat = msat + 1;
                else
                    break
                end
            end
            for j=1:nR % check row inequilities
                if Rineq(j,:)*L(:,Ridx(j)) <= 0
                    nRsat = nRsat + 1;
                else
                    nRsat = 0;
                    break
                end
            end
            if msat == m && nRsat == nR
                count_sat = count_sat + 1;
            end
        end
        if (mod(count_total, 10000) == 0)
            disp(['Out of ' num2str(count_total) ' posterior draws, ' num2str(count_sat) ' satisfy all the restrictions.']);
        end
    end
    Results_RWZ = [count_total/10^6, count_sat]*10000/time_set
    %%% ------------------proposed method-------------%%%
	Rineq = Rineq_newalg; % setup the row inequalities
    time_cost = tic; count_sat = 0; count_total = 0; k = 1;
    while toc(time_cost) < time_set
        [store_alp,store_beta,store_Sig] = sample_ThetaSig(Y0,Y,p,prior_redu,nbatch);
        count_total = count_total + nbatch;
        [store_Btilde,store_Sigtilde] = getReducedForm(store_alp,store_beta,store_Sig);
        for isim = 1:nbatch
            Sigtilde = squeeze(store_Sigtilde(isim,:,:));
            [Q,R] = QR(randn(n,n)); L0 = chol(Sigtilde,'lower'); L = L0*Q; satTab = zeros(m,n);
            for i=1:m  % check sign restrictions & row inequilities
                idx = find(S(:,i)==-1 | S(:,i)==1); nidx = length(idx); signL = sign(L(idx,:));
                for j=1:n
                    if (sum(signL(:,j) == S(idx,i)) == nidx) && ...
                            (sum((Rineq(i,:))*L(:,j) <= 0) == k) 
                        satTab(i,j) = 1;
                    elseif (sum(signL(:,j) == -S(idx,i)) == nidx) && ...
                            (sum((Rineq(i,:))*(-L(:,j)) <= 0) ==k) 
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
                reorder(m+1:end) = setdiff(1:n,reorder); L = L(:,reorder); % randomly permute and switch the signs of the last n-m columns
                tmpL2 = L(:,m+1:end); tmpL2 = tmpL2(:,randperm(n-m))*diag(2*(rand(n-m,1)>.5) - 1); L(:,m+1:end) = tmpL2;
                count_sat = count_sat + 1;
            end
        end
    end
    Results_proposed = [count_total/10^6, count_sat]*10000/time_set

    Tab3_PanelB(:,i_case) = [Results_RWZ, Results_proposed]
end
Tab3_PanelB = [Tab3_PanelB(:,1:3); Tab3_PanelB(:,4:6)]

%%% ---------------------------------- Panel C ----------------------------------%%%
Tab3_PanelC = zeros(4,6);
horizon = 36; dynamic_shock = [1,4,4,1,1,8]; K = 1; %pick the column of the dynamic shock and set the number of horizon K = 1
for i_case = 1:6
    eval(['load', ' SignRestr', num2str(i_case), '.mat']); eval(['load', ' DataT200_', num2str(i_case), '.mat']); % load data and sign restriction matrix
    eval(['load', ' RankRestr', num2str(i_case), '.mat']); % load ranking restriction matrix

    [T, n] = size(data); var_id = 1:n;
    Y0 = data(1:8,var_id); Y = data(9:end,var_id);
    [T,n] = size(Y); tmpY = [Y0(end-p+1:end,:); Y];
    Z = zeros(T,n*p);
    for ii=1:p
        Z(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
    end
    Z = [ones(T,1) Z]; kappa = [1, 1, 1, 100]; sig2 = get_resid_var(Y0,Y);
    prior_stru = prior_ACP_stru(n,p,kappa,sig2,[]); prior_redu = prior_ACP_redu(n,p,kappa,sig2,[]);

    S = Restr; m = size(S,2); % setup the sign restrictions
    %%% ------------------RWZ method-------------%%%
    Rineq = Rineq_RWZ; nR = length(Ridx); % setup the row inequalities	
    time_cost = tic; count_sat = 0; count_total = 0;
    while toc(time_cost) < time_set
        [store_alp,store_beta,store_Sig] = sample_ThetaSig(Y0,Y,p,prior_redu,nbatch);
        count_total = count_total + nbatch;
        [store_Btilde,store_Sigtilde] = getReducedForm(store_alp,store_beta,store_Sig);
        for isim = 1:nbatch
            Sigtilde = squeeze(store_Sigtilde(isim,:,:)); msat = 0; nRsat = 0;
            [Q,R] = QR(randn(n,n)); L0 = chol(Sigtilde,'lower'); L = L0*Q;
            for i=1:m  % check sign restrictions
                idx = find(S(:,i)==-1 | S(:,i)==1); nidx = length(idx); signL = sign(L(idx,:));
                if (sum(signL(:,i) == S(idx,i)) == nidx)
                    msat = msat + 1;
                elseif (sum(signL(:,i) == -S(idx,i)) == nidx)
                    L(:,i) = -L(:,i);
                    msat = msat + 1;
                else
                    break
                end
            end
            for j=1:nR % check row inequilities
                if Rineq(j,:)*L(:,Ridx(j)) <= 0
                    nRsat = nRsat + 1;
                else
                    nRsat = 0;
                    break
                end
            end
			Btilde = reshape(store_Btilde(isim,:),n*p+1,n); response = IRredu(Btilde(2:end,:),L,horizon,m); % includes IR at impact
            for i=dynamic_shock(i_case) % check dynamic restriction
                idx = find(S(:,i)==-1 | S(:,i)==1); nidx = length(idx); signL = sign(L(idx,:));
                sign_sat = sign(squeeze(response(idx,i,1:K+1))) == repmat(S(idx,i),1,K+1);
                if sum(sum(sign_sat)) == nidx*(K+1)
                    msat = msat + 1;
                else
                    break
                end
            end

            if msat == m && nRsat == nR
                count_sat = count_sat + 1;
            end
        end
        if (mod(count_total, 10000) == 0)
            disp(['Out of ' num2str(count_total) ' posterior draws, ' num2str(count_sat) ' satisfy all the restrictions.']);
        end
    end
    Results_RWZ = [count_total/10^6, count_sat]*10000/time_set
    %%% ------------------proposed method-------------%%%
	Rineq = Rineq_newalg; % setup the row inequalities
    time_cost = tic; count_sat = 0; count_total = 0; k = 1;
    while toc(time_cost) < time_set
        [store_alp,store_beta,store_Sig] = sample_ThetaSig(Y0,Y,p,prior_redu,nbatch);
        count_total = count_total + nbatch;
        [store_Btilde,store_Sigtilde] = getReducedForm(store_alp,store_beta,store_Sig);
        for isim = 1:nbatch
            Sigtilde = squeeze(store_Sigtilde(isim,:,:));
            [Q,R] = QR(randn(n,n)); L0 = chol(Sigtilde,'lower'); L = L0*Q; satTab = zeros(m,n);
            for i=1:m  % check sign restrictions & row inequilities
                idx = find(S(:,i)==-1 | S(:,i)==1); nidx = length(idx); signL = sign(L(idx,:));
                for j=1:n
                    if (sum(signL(:,j) == S(idx,i)) == nidx) && ...
                            (sum((Rineq(i,:))*L(:,j) <= 0) == k) 
                        satTab(i,j) = 1;
                    elseif (sum(signL(:,j) == -S(idx,i)) == nidx) && ...
                            (sum((Rineq(i,:))*(-L(:,j)) <= 0) ==k) 
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
                reorder(m+1:end) = setdiff(1:n,reorder); L = L(:,reorder); % randomly permute and switch the signs of the last n-m columns
                tmpL2 = L(:,m+1:end); tmpL2 = tmpL2(:,randperm(n-m))*diag(2*(rand(n-m,1)>.5) - 1); L(:,m+1:end) = tmpL2;
				Btilde = reshape(store_Btilde(isim,:),n*p+1,n); response = IRredu(Btilde(2:end,:),L,horizon,m);
                for i=dynamic_shock(i_case) % check dynamic sign restrictions (sign restrictions at impact are satisfied)
                    idx = find(S(:,i)==-1 | S(:,i)==1); nidx = length(idx);
                    sign_sat = sign(squeeze(response(idx,i,2:K+1))) == repmat(S(idx,i),1,K);
                    if sum(sum(sign_sat)) == nidx*K
                        msat = msat + 1;
                    else
                        break
                    end
                end
                if msat == 1
                    count_sat = count_sat + 1;
                end
            end
        end
    end
    Results_proposed = [count_total/10^6, count_sat]*10000/time_set

    Tab3_PanelC(:,i_case) = [Results_RWZ, Results_proposed]
end
Tab3_PanelC = [Tab3_PanelC(:,1:3); Tab3_PanelC(:,4:6)]

%%% ---------------------------------- Output for the whole table ----------------------------------%%%
Tab3 = [Tab3_PanelA; Tab3_PanelB; Tab3_PanelC]
csvwrite('Table3.csv',Tab3)