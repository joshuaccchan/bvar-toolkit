function qDraws = drawqs(restr,phi,opt)
% Function draws values of q satisfying sign and zero restrictions.
% Inputs:
% - restr: structure containing information about restrictions
% - phi: structure containing reduced-form VAR parameters
% - opt: structure containing model information and options

signRestr = restr.signRestr;
eqRestr = restr.eqRestr;
Sigmatr = phi.Sigmatr;
Sigmatrinv = phi.Sigmatrinv;
vma = phi.vma;

n = size(Sigmatrinv,1); % Number of variables in VAR
s = size(signRestr,1); % Number of sign restrictions
f = size(eqRestr,1); % Number of zero restrictions

%% Construct matrix representing zero restrictions.
% Restrictions are represented as F(phi)*q = 0.

F = zeros(f,n);

for ii = 1:f

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

S = zeros(s,n);

for ii = 1:s % For each restriction
    
    if signRestr(ii,5) == 1 % Sign restriction on impulse response
    
        S(ii,:) = vma(signRestr(ii,1),:,signRestr(ii,3)+1)*signRestr(ii,4);
    
    elseif signRestr(ii,5) == 2 % Sign restriction on A0
        
        S(ii,:) = Sigmatrinv(:,signRestr(ii,2))'*signRestr(ii,4);
        
    end
    
end

%% Draw values of q satisfying sign and zero restrictions.

qDraws = zeros(n,opt.L);

for kk = 1:opt.L

    flag = 0;

    while flag == 0
    
    % Draw nx1 vector of independent standard normal random variables.
    z = randn(n,1);
    
    % Compute residual from linear projection of z on F'.
    [q,r] = qr(F',0);
    q = z-F'*(r \ (q'*z));

    q = sign(Sigmatrinv(:,1)'*q)*(q./norm(q)); % Impose normalisations    

    % Check if sign restrictions satisfied given draw of q. If so, save
    % value of q and terminate while loop.
    flag = all(S*q >= 0);
    if flag == 1
        qDraws(:,kk) = q;
    end
    
    end
        
end

end