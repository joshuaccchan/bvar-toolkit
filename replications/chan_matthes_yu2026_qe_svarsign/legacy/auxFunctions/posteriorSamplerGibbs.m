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

%% Conduct posterior inference on impulse responses.
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
irfDraws = zeros(opt.phiDraws,opt.H+1,n);

nonEmptyDraw = 0; % Counter for no. of draws with nonempty identified set

tic
while nonEmptyDraw <= opt.phiDraws

    % Posterior sampler given normal-inverse-Wishart prior.
    phi.Sigma = iwishrnd(PhiTilde,nuTilde);
    phi.Sigmatr = chol(phi.Sigma,'lower');
    phi.Sigmatrinv = phi.Sigmatr\eye(n);
    phi.B = kron(phi.Sigmatr,cholOmegaTilde)*randn(m*n,1) + ...
        reshape(PsiTilde,n*m,1);
    
    % Generate coefficients in orthogonal reduced-form VMA representation.
    [phi.vma,~] = genVMA(phi,opt);
    
    % Check if identified set is empty and, if not, obtain initial value
    % for Gibbs sampler.
    [c0,K,restr.F,restr.S,restr.Sbar,empty] = chebyCheck(restr,phi,opt);
    
    if empty == 0 % If identified set is nonempty
    
        nonEmptyDraw = nonEmptyDraw + 1;
        
        % Store parameters.
        B(:,nonEmptyDraw) = phi.B;
        Sigmatr(:,:,nonEmptyDraw) = phi.Sigmatr;
        Sigmatrinv(:,:,nonEmptyDraw) = phi.Sigmatrinv;
        vma(:,:,:,nonEmptyDraw) = phi.vma;
        
        % Use Gibbs sampler to draw q from a uniform distribution over
        % the subspace satisfying the identifying restrictions.
        q0 = drawqGibbs(restr.Sbar,size(restr.F,1),K,c0,opt);
        
        % Compute impulse responses.
        for hh = 1:opt.H+1 % For each horizon
    
            % Compute vector of impulse responses.
            irfDraws(nonEmptyDraw,hh,:) = phi.vma(:,:,hh)*q0;

        end
    
    end

end
runTime = toc