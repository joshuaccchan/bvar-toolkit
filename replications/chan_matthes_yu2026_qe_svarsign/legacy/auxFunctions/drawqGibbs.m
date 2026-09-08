function q = drawqGibbs(Sbar,r,K,c0,opt)
% Function draws q from a uniform distribution over the subspace 
% satisfying the identifying restrictions using a Gibbs sampler.
% Inputs:
% - Sbar: matrix containing sign restrictions projected into linear
%   subspace satisfying zero restrictions (after change of basis).
% - r: number of zero restrictions
% - K: change of basis matrix
% - c0: initial value satisfying restrictions in transformed basis
% - opt: structure containing options

n = length(c0); % Dimension of linear subspace satisfying zero restrictions
L = opt.L; % Target number of draws
burn = opt.burn; % Number of initial draws to drop
f = opt.f; % Thinning parameter - keep every f-th draw.

z = zeros(n,burn + L*f);
z(:,1) = c0; % Initial value

% Obtain draws from a truncated multivariate normal distribution using a
% Gibbs sampler.

for kk = 2:(burn + L*f) % For each draw of z
    
    for jj = 1:n % For each element of z
        
        % Compute lower and upper bound of univariate truncated normal
        % distribution of z_{1,jj}.
        W = -(Sbar(:,1:jj-1)*z(1:jj-1,kk) + Sbar(:,jj+1:end)*z(jj+1:end,kk-1))...
            ./Sbar(:,jj);
        lb = max(W(Sbar(:,jj) > 0));
        if isempty(lb)
            lb = -Inf;
        end
        ub = min(W(Sbar(:,jj) < 0));
        if isempty(ub)
            ub = Inf;
        end      

        % Draw q_{1,jj} from truncated standard normal distribution using
        % inverse CDF.
        Fub = normcdf(ub);
        Flb = normcdf(lb);
        z(jj,kk) = norminv(rand*(Fub - Flb) + Flb,0,1);
              
    end
    
end

z = z(:,(burn+1):end); % Drop first burn draws
z = z(:,f:f:(L*f)); % Keep every f-th draw
% Project draws onto unit sphere.
q = z./vecnorm(z);
% Lift draws into n dimensions.
q = [q', zeros(L,r)]';
% Apply change of basis.
q = K*q;

end