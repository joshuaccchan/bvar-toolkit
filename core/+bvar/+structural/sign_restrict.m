% bvar.structural.sign_restrict - test one candidate rotation of the impact
% matrix against a set of sign restrictions and row inequalities, flipping the
% sign of any column that satisfies them in reverse. This is the accept-reject
% step of Rubio-Ramirez, Waggoner and Zha (2010), applied to rotations that
% bvar.structural.qr_sign draws uniformly from the orthogonal group.
%
%   [ok,L] = bvar.structural.sign_restrict(L, S, Rineq, Ridx)
%
%   L     : n x n candidate impact matrix, chol(Sigtilde,'lower')*Q for a random
%           rotation Q (see bvar.structural.qr_sign)
%   S     : n x m sign restrictions, one column per shock. +1 and -1 restrict the
%           sign of that response on impact; NaN leaves it free
%   Rineq : nR x n rows, each a linear combination required to be NEGATIVE
%   Ridx  : nR x 1 column of L that each row of Rineq applies to
%   ok    : true only if every shock satisfies its sign column AND every row
%           inequality holds
%   L     : the candidate with columns sign-flipped where that was what made the
%           restriction hold; unchanged when ok is false
%
% THE ROW INEQUALITIES ARE TESTED STRICTLY: every row must give
% Rineq(j,:)*L(:,Ridx(j)) < 0. A row of zeros therefore fails and rejects every
% candidate, producing an empty identified set with no error raised. To impose
% no ranking restrictions, pass an empty Ridx and a 0 x n Rineq, NOT a zero row.
% bvar.structural.sign_assign tests the same quantity as <= 0, where a zero row
% is harmless, so the two differ on exactly this input.
%
% A sign restriction identifies a shock only up to sign, so a column that
% violates S may satisfy it after negation, and the negated column is then the
% economically meaningful one. The check exits at the first shock that satisfies
% neither, which is why acceptance rates in the caller's rejection loop can be
% very low without any single evaluation being expensive.
%
% See:
% Rubio-Ramirez, J.F., Waggoner, D.F. and Zha, T. (2010). Structural Vector
% Autoregressions: Theory of Identification and Algorithms for Inference,
% Review of Economic Studies, 77(2): 665-696.
% Chan, J.C.C. (2022). Asymmetric Conjugate Priors for Large Bayesian VARs,
% Quantitative Economics, 13(3): 1145-1169.

function [ok,L] = sign_restrict(L,S,Rineq,Ridx)
m = size(S,2);
nR = numel(Ridx);
msat = 0;   % counter for the # of shocks that satisfies the sign restrictions
nRsat = 0;  % counter for the # of satisfied row inequalities

for i=1:m  % check sign restrictions
    idx = find(S(:,i)==-1 | S(:,i)==1);
    nidx = length(idx);
    signL = sign(L(idx,:));
        % check if the i-th column satisfies the sign restrictions
    if (sum(signL(:,i) == S(idx,i)) == nidx)
        msat = msat + 1;
        % or if the negative of the i-th column satisfies the sign restrictions
    elseif (sum(signL(:,i) == -S(idx,i)) == nidx)
        L(:,i) = -L(:,i); % change the sign of the i-th column
        msat = msat + 1;
    else
        break
    end
end
for j=1:nR % check row inequalities
    if Rineq(j,:)*L(:,Ridx(j)) < 0
        nRsat = nRsat + 1;
    else
        nRsat = 0;
        break
    end
end
ok = (msat == m && nRsat == nR);
end
