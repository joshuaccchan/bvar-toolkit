% This is the main script for replicating the application in Uhlig (2005) using three algorithms:
% algom = 1: RWZ (standard accept-reject algorithm);
% algom = 2: proposed algorithm;
% algom = 3: Read(2022) algorithm
% NOTICE: In compliance with licensing agreements, the commodity prices 
% presented have been perturbed with noise; the empirical results are virtually unaffected.

p = 12;        % if p > 12, need to change Y0 and Y below
model = 2;     % 1: VAR with ACP; 2: VAR with NCP
nsim = 5000;   % # of posterior draws (that satisfy all the restrictions)
K = 5;         % # of periods sign restrictions are imposed
addpath('./data');
addpath('./utility');
if model ==1
    nbatch = 10000; % # of posterior draws sampled in a batch
elseif model == 2
    nbatch = 500;
end

horizon = 61;   % # of horizons for impulse responses (including on impact)
% load data
data = readmatrix('./data/Uhlig_monthly.csv','Range','B602:G1069');
data(:,1:5) = log(data(:,1:5));
% 6 variables + 1 shock (monetary)
var_id = 1:6; %  GDP, inflation, commodity prices, nonborrowed reserve, total reserve, FFR
idx_ns = var_id;  % index for variables in levels
Y0 = data(1:12,var_id);  % save the first 12 obs as the initial conditions
Y = data(13:end,var_id);
[T,n] = size(Y);
tmpY = [Y0(end-p+1:end,:); Y];
Z = zeros(T,n*p);
for ii=1:p
    Z(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
end
Z = [ones(T,1) Z];
sig2 = get_resid_var(Y0,Y);
kappa = [.04,.0016,1,100];
if model == 1
    % find the optimal kappa values
    [ml_opt,kappa] = get_OptKappa(Y0,Y,Z,p,[.04,.0016],'redu',idx_ns);
    prior_stru = prior_ACP_stru(n,p,kappa,sig2,idx_ns);
    prior_redu = prior_ACP_redu(n,p,kappa,sig2,idx_ns);
elseif model == 2
    [prior.B0,prior.VB,prior.nu0,prior.S0] = prior_NCP(p,kappa(1),kappa(4),Y0,Y);
    prior.B0 = [zeros(1,n);sparse(idx_ns,idx_ns,ones(1,length(idx_ns)),n,n);zeros((p-1)*n,n)];  % random walk prior mean
end
% setup the sign restrictions and row inequalities
monetary = [0,-1,-1,-1,0,1]';
S = monetary;
m = size(S,2); % # of shocks

start_time = clock;
count_sat = 0; % counter for # draws that satisfy all the conditions
count_total = 0; % counter for total # draws
checksigns = 0;
disp(['Computing impulse responses from a ' num2str(n) '-variable VAR']);
disp('    to an one-standard-deviation monetary policy shock...');
store_response = zeros(nsim,n,m,horizon);
while count_sat < nsim
    switch model
        case 1
            % sample nbatch draws from the posterior
            [store_alp,store_beta,store_Sig] = sample_ThetaSig(Y0,Y,p,prior_redu,nbatch);
            % obtain the reduced-form parameters
            [store_Btilde,store_Sigtilde] = getReducedForm(store_alp,store_beta,store_Sig);
        case 2
            [store_Btilde,store_Sigtilde] = sample_BSig_NCP(Y0,Y,p,prior,nbatch);
    end
    count_total = count_total + nbatch;

    for isim = 1:nbatch  % go through the nbatch posterior draws to find those
        Sigtilde = squeeze(store_Sigtilde(isim,:,:));
        msat = 0;  % counter for the # of shocks that satisfies the sign restrictions

        [Q,R] = QR(randn(n,n));
        L0 = chol(Sigtilde,'lower');
        L = L0*Q;

        switch algom
            case 1 % use standard accept-reject
                Btilde = reshape(store_Btilde(isim,:),n*p+1,n);
                response = IRredu(Btilde(2:end,:),L,horizon,m); % includes IR at impact
                for i=1:m
                    idx = find(S(:,i)==-1 | S(:,i)==1);
                    nidx = length(idx);
                    signL = sign(L(idx,:));
                    % check sign restrictions
                    sign_sat = sign(squeeze(response(idx,i,1:K+1))) == repmat(S(idx,i),1,K+1);
                    if sum(sum(sign_sat)) == nidx*(K+1)
                        msat = msat + 1;
                    else
                        break
                    end
                end
                if msat == m
                    count_sat = count_sat + 1;
                    store_response(count_sat,:,:,:) = response;
                end
            case 2 % use proposed algorithm
                satTab = zeros(m,n); % (i,j) = 1 if the j-th column of L satisfies all restrictions for the i-th shock
                % (i,j) = -1 if the negative of j-th column of L satisfies all restrictions
                for i=1:m  % check sign restrictions & row inequilities
                    idx = find(S(:,i)==-1 | S(:,i)==1);
                    nidx = length(idx);
                    signL = sign(L(idx,:));
                    for j=1:n
                        if (sum(signL(:,j) == S(idx,i)) == nidx)
                            satTab(i,j) = 1;
                        elseif (sum(signL(:,j) == -S(idx,i)) == nidx)
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
                    checksigns = checksigns + sum(sum(S.*L<0));  % count # of violations
                    % randomly permute and switch the signs of the last n-m columns
                    tmpL2 = L(:,m+1:end);
                    tmpL2 = tmpL2(:,randperm(n-m))*diag(2*(rand(n-m,1)>.5) - 1);
                    L(:,m+1:end) = tmpL2;
                    Btilde = reshape(store_Btilde(isim,:),n*p+1,n);
                    response = IRredu(Btilde(2:end,:),L,horizon,m);

                    % check dynamic sign restrictions (sign restrictions at impact are satisfied)
                    for i=1:m
                        idx = find(S(:,i)==-1 | S(:,i)==1);
                        nidx = length(idx);
                        sign_sat = sign(squeeze(response(idx,i,2:K+1))) == repmat(S(idx,i),1,K);
                        if sum(sum(sign_sat)) == nidx*K
                            msat = msat + 1;
                        else
                            break
                        end
                    end
                    if msat == m
                        count_sat = count_sat + 1;
                        store_response(count_sat,:,:,:) = response;
                    end
                end
            case 3 % use Read(2022) algorithm
                %  GDP:0, inflation:-1, commodity prices:-1, nonborrowed reserve:-1, total reserve:0, FFR:+1
                %% Baseline set of sign restrictions.
                % Each row of signRestr contains a vector (i,j,h,s,t) representing a sign
                % restriction, where t is the type of restriction:
                % t = 1: the response of the ith variable with respect to the jth shock at
                % the hth horizon is nonnegative (s = 1) or nonpositive (s = -1).
                % t = 2: the (ij)th element of A0 is nonnegative (s = 1) or nonpositive
                % (s = -1).
                H = 5;
                S1 = [2 1 0 -1 1]; % Impact response of inflation to MPS is nonpositive
                S1 = repmat(S1,[H+1,1]);
                S1(:,3) = 0:H;
                S2 = [3 1 0 -1 1]; % Impact response of commodity prices to MPS is nonpositive
                S2 = repmat(S2,[H+1,1]);
                S2(:,3) = 0:H;
                S3 = [4 1 0 -1 1]; % Impact response of nonborrowed reserve to MPS is nonpositive
                S3 = repmat(S3,[H+1,1]);
                S3(:,3) = 0:H;
                S4 = [6 1 0 1 1]; % Impact response of fed funds rate to MPS is nonpositive
                S4 = repmat(S4,[H+1,1]);
                S4(:,3) = 0:H;
                restr.signRestr = [S1; S2; S3; S4];
                restr.eqRestr = [];

                %% Options.
                opt.p = 12; % No. of lags in VAR
                opt.const = 1; % const = 1 if constant in VAR, = 0 otherwise
                opt.nExog = opt.const;
                opt.cumIR = []; % Indices of variables for which cumulative impulse response should be used
                opt.H = 60; % Terminal horizon for impulse responses
                opt.burn = 3; % Drop first burn draws from Gibbs sampler
                opt.L = 1;
                opt.f = 1;
                opt.phiDraws = 1000; % No. of draws from posterior of phi

                addpath('auxFunctions');
                % use Read's notation
                phi.Sigma = Sigtilde;
                phi.Sigmatr = L0;
                phi.Sigmatrinv = phi.Sigmatr\eye(n);
                Btilde = reshape(store_Btilde(isim,:),n*p+1,n);
                phi.B = store_Btilde(isim,:)';

                % Generate coefficients in orthogonal reduced-form VMA representation.
                [phi.vma,~] = genVMA(phi,opt);

                % Use Algorithm 1 to determine whether identified set is nonempty and,
                % if so, obtain a value of q satisfying the sign and zero restrictions
                % in the transformed basis (i.e., in an (n-r)-dimensional subspace).
                [c0,K,restr.F,restr.S,restr.Sbar,empty] = chebyCheck_nonorm(restr,phi,opt);

                if empty == 0 % If identified set is nonempty

                    %disp(['The draw from ' num2str(count_sat) ' is nonempty.']);
                    q0 = drawqGibbs(restr.Sbar,size(restr.F,1),K,c0,opt);
                    L = L0*[q0, zeros(n,n-1)]; %L(:,1)'
                    response = IRredu(Btilde(2:end,:),L,horizon,m);
                    count_sat = count_sat + 1;
                    store_response(count_sat,:,:,:) = response;

                end
        end
    end

    if (mod(count_total, 1000) == 0)
        disp(['Out of ' num2str(count_total) ' posterior draws, ' ...
            num2str(count_sat) ' satisfy all the restrictions.']);
    end
end
disp([num2str(nsim) ' posterior draws that satisfy all the restrictions are obtained.']);
disp(['The simulation took ' num2str(etime(clock,start_time)/60) ' minutes.']);
