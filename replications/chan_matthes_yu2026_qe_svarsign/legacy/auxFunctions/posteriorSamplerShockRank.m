tic

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
ni = length(opt.ivar); % Number of variables of interest
m = size(XX,2); % Number of parameters in each equation for y_t
opt.nExog = opt.const + size(exog,2); % Number of exogenous variables
T = length(YY); % Number of observations used in estimating VAR

% Adjust indices representing narrative restrictions to account for
% losing the first p observations.
if ~isempty(restr.shockRank)
    restr.shockRank(:,1) = restr.shockRank(:,1) - opt.p;
end

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

irfDraws = zeros(opt.Kdraws,opt.H+1,length(opt.ivar));
irfMin = irfDraws;
irfMax = irfDraws;

phiDraw = 0; % Counter for no. of draws from posterior of phi
draw = 0; % Counter for no. of draws with non-empty IS

tic

while draw < opt.Kdraws

    % Posterior sampler given normal-inverse-Wishart prior.
    phiDraw = phiDraw + 1;
    phi.Sigma = iwishrnd(PhiTilde,nuTilde);
    phi.Sigmatr = chol(phi.Sigma,'lower');
    phi.Sigmatrinv = phi.Sigmatr\eye(n);
    phi.B = kron(phi.Sigmatr,cholOmegaTilde)*randn(m*n,1) + ...
        reshape(PsiTilde,n*m,1);
    
    % Generate coefficients in orthogonal reduced-form VMA representation.
    [phi.vma,~] = genVMA(phi,opt);
    
    % Compute reduced-form VAR innovations.
    restr.U = (YY - XX*reshape(phi.B,m,n))';
        
    % Use Algorithm 1 to determine whether identified set is nonempty and,
    % if so, obtain a value of q satisfying the sign and zero restrictions 
    % in the transformed basis (i.e., in an (n-r)-dimensional subspace).
    [c0,K,restr.F,restr.S,restr.Sbar,Qempty(phiDraw)] = ...
        chebyCheck_shockRank(restr,phi,opt);
               
    if Qempty(phiDraw) == 0 % If (conditional) identified set is non-empty
        
            draw = draw + 1;
            
            % Use Gibbs sampler to draw q from a uniform distribution over
            % the subspace satisfying the identifying restrictions.
            q0 = drawqGibbs(restr.Sbar,size(restr.F,1),K,c0,opt);
            
            for hh = 1:opt.H+1 % For each horizon
    
                % Compute vector of impulse responses using final draw of 
                % q from Gibbs sampler.
                irfDraws(draw,hh,:) = phi.vma(opt.ivar,:,hh)*q0(:,opt.L);

            end
            
            % Transform value of q from Algorithm 4.1 into original basis.
            c0 = K*[c0', zeros(1,size(restr.F,1))]';
                
            % Compute identified sets for impulse responses of interest.
           [irfMin(draw,:,:),irfMax(draw,:,:)] = ...
               numericalBounds(restr,phi,c0,opt,optimOptions);  
                      
           if mod(opt.Kdraws-draw,opt.dispIter) == 0
              
               fprintf('\n%d draws with non-empty identified set remaining...',...
               opt.Kdraws-draw);
               
           end

    end

end

runTime = toc

% Compute cumulative impulse response where necessary.
if ~isempty(opt.cumIR)
  
    irfDraws(:,:,opt.cumIR) = cumsum(irfDraws(:,:,opt.cumIR),2);

end

% Compute posterior mean under single prior.
irfMean = permute(mean(irfDraws,1),[2 3 1]);

% Compute highest posterior density intervals under single prior.
[irflb,irfub] = highestPosteriorDensity(irfDraws,opt);

% Compute set of posterior means.
irfMeanlb = permute(mean(irfMin,1),[2 3 1]);
irfMeanub = permute(mean(irfMax,1),[2 3 1]);

% Compute robust credible intervals.
[irfCredlb,irfCredub] = credibleRegion(irfMin,irfMax,opt);

% Compute posterior plausibility of the identifying restrictions.
postPlaus = draw/phiDraw;
fprintf('\nPosterior plausibility of identifying restrictions: %0.4g.\n',...
    postPlaus)
