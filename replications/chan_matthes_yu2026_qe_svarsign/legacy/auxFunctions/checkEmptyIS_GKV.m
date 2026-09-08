function [empty,F,S] = checkEmptyIS_GKV(restr,phi,opt)
% Function determines whether the identified set is empty using the
% approach in Giacomini, Kitagawa and Volpicella (2020).
% Specifically, consider all possible combinations of n-r-1 active (or 
% binding) sign restrictions and check whether there exists a vector
% satisfying these constraints and any inactive constraints.
% Inputs:
% - restr: structure containing information about restrictions
% - phi: structure containing reduced-form VAR parameters
% - opt: structure containing model information and options

signRestr = restr.signRestr;
eqRestr = restr.eqRestr;
Sigmatr = phi.Sigmatr;
Sigmatrinv = phi.Sigmatrinv;
vma = phi.vma;

n = size(Sigmatr,1); % Number of variables in VAR
s = size(signRestr,1); % Number of sign restrictions (excluding normalisation)
r = size(eqRestr,1); % Number of zero restrictions

%% Construct matrix representing zero restrictions.
% Restrictions are represented as F(phi)*q = 0.

F = zeros(r,n);

for ii = 1:r

    if eqRestr(ii,3) == 1 % Restriction on A0

        F(ii,:) = Sigmatrinv(:,eqRestr(ii,2))';

    elseif eqRestr(ii,3) == 2 % Restriction on A0^(-1)

        F(ii,:) = Sigmatr(eqRestr(ii,1),:);

    elseif eqRestr(ii,3) == 3 % Restriction on LRCIR

        % Reshape VAR coefficients into matrix.
        B = reshape(phi.B,size(phi.B,1)/n,n);
        B = B(1:end-opt.nExog,:); % Drop coefficients on exogenous variables
        B = reshape(B',[n,n,opt.p]); % Reshape into array

        % Compute matrix of reduced-form LRCIRs post-multiplied by Sigma_tr.
        lrcir = ((eye(n)-sum(B,3))\eye(N))*Sigmatr;
    
        % Extract relevant row.
        F(ii,:) = lrcir(eqRestr(ii,1),:);

    end

end

%% Construct matrix representing sign restrictions.
% Restrictions are represented as S(phi)*q >= 0.

S = zeros(s+1,n);

for ii = 1:s % For each restriction
    
    if signRestr(ii,5) == 1 % Sign restriction on impulse response
    
        S(ii,:) = vma(signRestr(ii,1),:,signRestr(ii,3)+1)*signRestr(ii,4);
    
    elseif signRestr(ii,5) == 2 % Sign restriction on A0
        
        S(ii,:) = Sigmatrinv(:,signRestr(ii,2))'*signRestr(ii,4);
        
    end
    
end

S(end,:) = Sigmatrinv(:,1)'; % Add sign normalisation

s = s + 1; % Redefine number of sign restrictions to include normalisation

%% If s + r <= n, identified set is never empty.
% This bit is commented out so that comparison between algorithms is fair.
% Uncomment if using in applications.
% if s + r <= n
%     empty = 0;
%     return
% end

%% Check for vector satisfying restrictions.

% Number of active sign restrictions to consider
sTilde = n-r-1;

% Generate matrix containing indices representing all possible 
% combinations of sTilde of the columns of S.
sComb = nchoosek(1:s,sTilde);
% Each row of sComb contains n-r-1 indices representing a combination of
% sign restrictions (indexing the rows of S).
rNum = size(sComb,1); % Number of possible combinations

for rr = 1:rNum % For each possible combination of active constraints

    % Construct matrix containing active constraints.
    Z = [F; S(sComb(rr,:),:)];
    % Construct matrix containing non-active sign restrictions.
    Sna = S(setdiff(1:end,sComb(rr,:)),:);
    % Compute an orthonormal basis for null space of Z, which is a
    % one-dimensional vector since rank(Z) = n-1.
    q = null(Z);
    % Check if q or -q satisfy the non-active sign restrictions.
    empty = ~(all(Sna*q >= 0) || all(-Sna*q >= 0));

    if empty == 0 % If set is nonempty, terminate function

        return

    end

end  
        
end