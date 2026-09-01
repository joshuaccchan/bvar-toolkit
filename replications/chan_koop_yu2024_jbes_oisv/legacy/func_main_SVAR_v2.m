% This is the main run file for estimating the large Bayesian VAR (OI, CS)
% with the Minnesota-type horseshoe prior
%
function [Mean_kappa,Mean_Impact,Median_kappa,Median_Impact,Sig_mean,Sig_median] = func_main_SVAR_v2(model,rev_option,nsim,burnin)
addpath('./utility')
%nsim = 10000;
%burnin = 1000;
p = 13; % if p > 8, need to change Y0 and Y below

% 1: Minnesota-type Horseshoe SVARSV_redu; 2: Minnesota-type Horseshoe CCCM_redu;
%model = 1; rev_option = 0; % rev_option-0: default; rev_option-1: reversed order;

% load data
data = load('FRED_MD_20vars.csv'); %4--IP, 6--unemployment, 12--PCE inflation, 13--FEDFUNDS
var_id1 = [4,6,12,13]; % n = 23
var_id2 = [1:3,5,7:11,14:20]; % first four variables: GDP growth, unemployment, PCE inflation, FEDFUNDS rate
var_id = [var_id1, var_id2];

if rev_option== 1
    var_id = flip(var_id);
end
Y0 = data(1:24,var_id);  % save the first 24 obs as the initial conditions
Y = data(25:end,var_id);

[T,n] = size(Y);
tmpY = [Y0(end-p+1:end,:); Y];
X = zeros(T,n*p);
for ii=1:p
    X(:,(ii-1)*n+1:ii*n) = tmpY(p-ii+1:end-ii,:);
end
X = [ones(T,1) X];

switch model
    case 1
        SVARSV_MH;
        Mean_kappa = kappa_mean; Median_kappa = median(store_kappa)';
        Mean_Impact = reshape(B0_mean,n,n)';
        B0_median = median(store_B0)'; Median_Impact = reshape(B0_median,n,n)';
        
        store_Sig = zeros(T,n,n);
        for isim = 1:nsim
            B0 = reshape(store_B0(isim,:),n,n)'; % stacked by rows
            h = squeeze(store_h(isim,:,:));
            store_Sig = store_Sig + construct_Sigt(h,B0);
        end
        Sig_mean = store_Sig/nsim; Sig_median = zeros(T,n,n);       
    case 2
        CS_MH;
        Mean_kappa = kappa_mean; Median_kappa = median(store_kappa)';
        A_id = nonzeros(tril(reshape(1:n^2,n,n),-1)'); A = eye(n); A(A_id) = alp_mean;
        Mean_Impact = A;
        A_id = nonzeros(tril(reshape(1:n^2,n,n),-1)'); A = eye(n); alp_median = median(store_alp)'; A(A_id) = alp_median;
        Median_Impact = A;
        
        store_Sig = zeros(T,n,n);
        for isim = 1:nsim
            h = squeeze(store_h(isim,:,:));
            alp = squeeze(store_alp(isim,:))';
            A = eye(n); A(A_id) = alp;
            store_Sig = store_Sig + construct_Sigt(h,A);
        end
        Sig_mean = store_Sig/nsim; Sig_median = zeros(T,n,n);        
end
end


function Sigt = construct_Sigt(h,B0)
[T,n] = size(h);
Sigt = zeros(T,n,n);
for t=1:T
    Sigt(t,:,:) = (B0\sparse(1:n,1:n,exp(h(t,:))))/(B0');
end
end