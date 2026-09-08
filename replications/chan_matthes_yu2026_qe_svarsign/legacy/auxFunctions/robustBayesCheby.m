%% Construct data for estimating VAR for y_t.
YY = data(opt.p+1:end,:); % y_t 
XX = lagmatrix(data,1:opt.p); % Matrix of regressors in VAR for y_t
XX = XX(opt.p+1:end,:); % Drop initial missing observations
if opt.const == 1 % Add constant to matrix of regressors
    XX = [XX ones(size(XX,1),1)];
end
% Add exogenous variables to matrix of regressors
XX = [XX exog(opt.p+1:end,:)]; 

n = size(YY,2); % Number of variables
m = size(XX,2); % Number of parameters in each equation for y_t
opt.nExog = opt.const + size(exog,2); % Number of exogenous variables
T = length(YY); % Number of observations used in estimating VAR

%% Conduct robust Bayesian inference on impulse responses.
% Assumes normal-inverse-Wishart prior.
nuBar = 0; % Prior degrees of freedom
PhiBar = zeros(n); % Prior scale
PsiBar = zeros(m,n); % Prior location of regression coefficients
% Inverse of prior variance of regression coefficients.
OmegaBarInv = zeros(m);

% Derive posterior parameters.
nuTilde = T + nuBar; % Degrees of freedom
OmegaTilde = (XX'*XX + OmegaBarInv)\eye(m); % Variance.
cholOmegaTilde = chol(OmegaTilde,'lower');
OmegaTildeInv  =  XX'*XX  + OmegaBarInv;
PsiTilde = OmegaTilde*(XX'*YY + OmegaBarInv*PsiBar); % Location
PhiTilde = YY'*YY + PhiBar + PsiBar'*OmegaBarInv*PsiBar - ...
    PsiTilde'*OmegaTildeInv*PsiTilde; % Scale
PhiTilde = (PhiTilde+PhiTilde')*0.5; % Ensure matrix symmetric

% Storage arrays.
B = zeros(n*m,opt.phiDraws);
Sigmatr = zeros(n,n,opt.phiDraws);
Sigmatrinv = zeros(n,n,opt.phiDraws);
vma = zeros(n,n,opt.H+1,opt.phiDraws);
irfMin = zeros(opt.phiDraws,opt.H+1,length(opt.ivar));
irfMax = zeros(opt.phiDraws,opt.H+1,length(opt.ivar));

nonEmptyDraw = 0; % Counter for no. of draws with nonempty identified set
phiDraw = 0;

tic
while nonEmptyDraw < opt.phiDraws

    % Posterior sampler given normal-inverse-Wishart prior.
    phiDraw = phiDraw+1;
    phi.Sigma = iwishrnd(PhiTilde,nuTilde);
    phi.Sigmatr = chol(phi.Sigma,'lower');
    phi.Sigmatrinv = phi.Sigmatr\eye(n);
    phi.B = kron(phi.Sigmatr,cholOmegaTilde)*randn(m*n,1) + ...
        reshape(PsiTilde,n*m,1);
    
    % Generate coefficients in orthogonal reduced-form VMA representation.
    [phi.vma,~] = genVMA(phi,opt);
    
    % Use Algorithm 1 to determine whether identified set is nonempty and,
    % if so, obtain a value of q satisfying the sign and zero restrictions 
    % in the transformed basis (i.e., in an (n-r)-dimensional subspace).
    [q0,K,restr.F,restr.S,restr.Sbar,empty] = chebyCheck(restr,phi,opt);
    
    if empty == 0 % If identified set is nonempty
    
        nonEmptyDraw = nonEmptyDraw + 1;
        
        % Transform value of q from Algorithm 4.1 into original basis.
        q0 = K*[q0', zeros(1,size(restr.F,1))]';
                
        % Compute identified sets for impulse responses of interest.
       [irfMin(nonEmptyDraw,:,:),irfMax(nonEmptyDraw,:,:)] = ...
           numericalBounds(restr,phi,q0,opt,optimOptions);  
  
    end

end
runTime = toc

% Compute set of posterior means.
irfMeanlb = permute(mean(irfMin,1),[2 3 1]);
irfMeanub = permute(mean(irfMax,1),[2 3 1]);

% Compute robust credible intervals.
[irfCredlb,irfCredub] = credibleRegion(irfMin,irfMax,opt);