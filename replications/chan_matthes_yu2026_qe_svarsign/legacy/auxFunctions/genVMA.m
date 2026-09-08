function [vma,nonStable] = genVMA(phi,opt)
% Given reduced-form VAR parameters, generate coefficients in orthogonal
% reduced-form vector moving average (VMA) representation and check 
% stability of the VAR representation.
% Inputs:
% - phi: structure containing reduced-form VAR parameters
% - opt: structure containing model information and options
% - nExog: number of exogenous variables.

Sigmatr = phi.Sigmatr;
B = phi.B;
n = size(Sigmatr,1); % Number of variables in VAR

%% Put VAR(p) in companion form (i.e. VAR(1) for (y_t,y_t-1,...,y_t-p+1)').
B_mat = reshape(B,size(B,1)/n,n); % Reshape coefficients into matrix
%B_mat = B_mat(1:end-opt.nExog,:); % Drop coefficients on exogenous variables
B_mat = B_mat(2:end,:); % Drop coefficients on exogenous variables %Yu modified 
B_c = zeros(n*opt.p);
B_c(1:n,:) = B_mat';   
B_c(n+1:end,1:end-n) = eye(n*opt.p-n);

%% Compute coefficients in vector moving average representation.
vma = zeros(n,n,opt.H+1); % Storage array 
CC = eye(n*opt.p);
vma(:,:,1) = Sigmatr;

for jj = 2:opt.H+1 % For each horizon

    CC = CC*B_c;
    % Extract upper-left nxn block of coefficient matrix.
    vma(:,:,jj) = CC(1:n,1:n)*Sigmatr;
    
end

%% Compute cumulative impulse responses for relevant variables.

if ~isempty(opt.cumIR) 
    
    vma(cumIR,:,:) = cumsum(vma(opt.cumIR,:,:),3);
    
end

%% Check stability of VAR.
% Compute eigenvalues of companion matrix.
E = eig(B_c);
% Check if any eigenvalues are greater than one in modulus.
nonStable = any(abs(E) >= 1);

end