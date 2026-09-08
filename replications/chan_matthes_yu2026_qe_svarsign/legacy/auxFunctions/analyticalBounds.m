function [rmin,rmax] = analyticalBounds(restr,phi,opt)
% Function uses analytical results from Gafarov, Meier and Montiel-Olea
% (2018) (GMM18) to find the upper and lower bounds of the impulse-response
% identified set.
% Inputs:
% - restr: structure representing restrictions
% - phi: structure containing reduced-form VAR parameters
% - opt: structure containing model information and options

F = restr.F;
S = restr.S;
vma = phi.vma;
Sigmatr = phi.Sigmatr;
Sigma = Sigmatr*Sigmatr';
Sigmatrinv = phi.Sigmatrinv;
ivar = opt.ivar;
H = opt.H;

n = size(Sigmatr,1); % Number of variables in VAR
r = size(F,1); % Number of zero restrictions
s = size(S,1); % Number of sign restrictions

% Convert restrictions to be consistent with notation from GMM18.
F = (F*Sigmatrinv)';
S = (S*Sigmatrinv)';

% Compute number of combinations of active sign restrictions.
kmax = min(s,n-r-1); % Maximum number of active sign restrictions
rNum = 0;

for kk = 0:kmax

    rNum = nchoosek(s,kk) + rNum;

end

% Storage
fMax = zeros(H+1,rNum,length(ivar));
fMin = zeros(H+1,rNum,length(ivar));

cbar = 10^8; % Set penalty term

%% No active sign restrictions (k=0).
Z = F; % Active restrictions (equality restrictions only)

M = eye(n)-Sigmatr'*Z*((Z'*Sigma*Z)\Z'*Sigmatr);

for hh = 1:H+1 % For each horizon

    for ii = 1:length(ivar) % For each variable of interest

        % Reduced-form IR of ith variable at hth horizon.
        CC = vma(ivar(ii),:,hh)*Sigmatrinv; 

        v = real(sqrt(CC*Sigmatr*M*Sigmatr'*CC')); % Value function           

        if v ~= 0 % If value function non-zero

            % Potential solutions.
            xPlus = Sigmatr*M*Sigmatr'*CC'./v;
            xMinus = -xPlus;

            fMaxPlus = v - 2*(1-all(S'*xPlus >= 0))*cbar;
            fMaxMinus = -v - 2*(1-all(S'*xMinus >= 0))*cbar;

            fMinPlus = v + 2*(1-all(S'*xPlus >= 0))*cbar;
            fMinMinus = -v + 2*(1-all(S'*xMinus >= 0))*cbar;

            fMax(hh,1,ii) = max(fMaxPlus,fMaxMinus);
            fMin(hh,1,ii) = min(fMinPlus,fMinMinus);          

        elseif v == 0
            
            % Check whether there exists a vector xstar!=0 satisfying the
            % active restrictions in Z and the nonactive sign restrictions.
            
            if size(Z,2) < n-1
                
                % Use Chebyshev criterion when number of active
                % restrictions is less than n-1. 
                empty = chebyCheck2(Z'*Sigmatr,S'*Sigmatr);
                
            elseif size(Z,2) == n-1
                
                % If number of restrictions is equal to n-1, the vector
                % satisfying the active restrictions is pinned down up to
                % sign. Compute unit-length vector satisfying active
                % restrictions.
                xstar = null(Z'*Sigmatr);
                empty = ~(all(S'*Sigmatr*xstar >= 0) | ...
                    all(-S'*Sigmatr*xstar >= 0));  
                
            end
            
            if empty == 0

                fMax(hh,1,ii) = 0;
                fMin(hh,1,ii) = 0;

            elseif empty == 1

                fMax(hh,1,ii) = -2*cbar;
                fMin(hh,1,ii) = 2*cbar;

            end

        end

    end

end


%% 0 < k <= kmax (for different numbers of active sign restrictions)

rCount = 1; % Counter for combinations of active restrictions

for kk = 1:kmax

    % Generate vector containing column indices representing all 
    % possible combinations of k of the columns of S.
    sComb = nchoosek(1:s,kk);

    for rr = 1:size(sComb,1) % For each set of active constraints

        rCount = rCount + 1;
        
        % Active restrictions.
        Z = [F S(:,sComb(rr,:))];
        % Non-active sign restrictions.
        Sna = S(:,setdiff(1:end,sComb(rr,:)));
        
        M = eye(n) - Sigmatr'*Z*((Z'*Sigma*Z)\Z'*Sigmatr);

        for hh = 1:H+1 % For each horizon   

            for ii = 1:length(ivar) % For each variable of interest
                
                % Reduced-form IR of ith variable at hth horizon.
                CC = vma(ivar(ii),:,hh)*Sigmatrinv;
                
                v = real(sqrt(CC*Sigmatr*M*Sigmatr'*CC')); % Value function
                
                if v~=0 % If value function non-zero
                
                    % Potential solutions.
                    xPlus = Sigmatr*M*Sigmatr'*CC'./v;
                    xMinus = -xPlus;

                    fMaxPlus = v - 2*(1-all(Sna'*xPlus >= 0))*cbar;
                    fMaxMinus = -v - 2*(1-all(Sna'*xMinus >= 0))*cbar;

                    fMinPlus = v + 2*(1-all(Sna'*xPlus >= 0))*cbar;
                    fMinMinus = -v + 2*(1-all(Sna'*xMinus >= 0))*cbar;

                    fMax(hh,rCount,ii) = max(fMaxPlus,fMaxMinus);
                    fMin(hh,rCount,ii) = min(fMinPlus,fMinMinus);
                    
                elseif v == 0
                    
                    % Check whether there exists a point xstar!=0 
                    % satisfying the active restrictions in Z and the 
                    % nonactive sign restrictions. 
                    
                    if size(Z,2) < n-1

                        % Use Chebyshev criterion when number of active
                        % restrictions is less than n-1. 
                        empty = chebyCheck2(Z'*Sigmatr,Sna'*Sigmatr);

                    elseif size(Z,2) == n-1

                        % If number of restrictions is equal to n-1, the 
                        % vector satisfying the active restrictions is 
                        % pinned down up to sign. Compute unit-length 
                        % vector satisfying active restrictions.
                        xstar = null(Z');
                        empty = ~(all(Sna'*Sigmatr*xstar >= 0) | ...
                            all(-Sna'*Sigmatr*xstar >= 0));  

                    end

                    if empty == 0

                        fMax(hh,1,ii) = 0;
                        fMin(hh,1,ii) = 0;

                    elseif empty == 1

                        fMax(hh,1,ii) = -2*cbar;
                        fMin(hh,1,ii) = 2*cbar;

                    end                    
               
                end

            end

        end

    end

end

% Compute bounds as max/min over fMax/fMin, respectively, for each
% horizon and variable.
rmax = permute(max(fMax,[],2),[1 3 2]);
rmin = permute(min(fMin,[],2),[1 3 2]);

end