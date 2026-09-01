function [b_0, B_0, v_0, S_0] = Min_Prior(data0, lag, ar_order, kepa_1, kepa_2, kepa_3)
global  n
    K = 1 + n * lag;
    Y = data0;   
    
    % To get the residual variances of univariate p-lag autoregressions, we
    % run the AR(p) model on each equation, ignoring the constant
    % and exogenous variables (if they have been specified for the original
    % VAR model)
    sigma_sq = zeros(n, 1); % vector to store residual variances
    I = zeros(n*K, 2);
    % AR regression
    for i=1:n
        T = size(Y(:,i),1);
        X = [ones(T-4,1) Y(4:end-1,i) Y(3:end-2,i) Y(2:end-3,i) Y(1:end-4,i)];
        ols = (X'*X)\(X'*Y(5:end,i));
        sigma_sq(i,1) = mean((Y(5:end,i) - X*ols).^2);
        % m = ar(Y(:,i), ar_order);
        % sigma_sq(i,1) = m.Report.Fit.MSE;
    end    
    % Get B_0 as a diagonal matrix
    V_i = zeros(n,K);
    V_i(:,1) = kepa_2 * ones(n,1); % the intercept terms
    I(1:n,2) = ones(n,1);
    V_i(:,2:end) = repmat(sigma_sq,1,K-1); % the other B depends on which variables 
    % Divide lag^2 to discount
    for j = 1:lag
        ind_rng = 2+(j-1)*n:n*j+1;
        I(n+(j-1)*n*n+1:n+j*n*n,1) = reshape(V_i(:,2+(j-1)*n:n*j+1)/j^2, n^2, 1); 
        V_i(:,ind_rng) = V_i(:,ind_rng) * kepa_1./j^2;  
    end    
 
   K_nq = commutation_matrix(n,K);
   I = K_nq * I;
   b_0.v = zeros(K*n, 1);   b_0.d=zeros(K*n,3);
   B_0.v = diag(reshape(V_i',K*n,1));   
   B_0.d=sparse([d_diag(I,K*n),zeros(K*n*K*n,1)]);
   v_0 = n + 3;
   S_0.v = kepa_3 * eye(n);   
   S_0.d=sparse([zeros(n*n,2),d_diag(ones(n,1),n)]);
  S_0=dmatrix(S_0,n,3);
B_0=dmatrix(B_0,K*n,3);
end
