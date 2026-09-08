% bvar.structural.sign_assign - accept a candidate impact matrix if EVERY shock
% has at least one column that satisfies its restrictions, then draw one such
% assignment at random. Sign restrictions and ranking restrictions are handled
% together, the first being the special case of the second with a zero weight.
% This is the fast alternative to bvar.structural.sign_restrict.
%
%   [ok,L] = bvar.structural.sign_assign(L, S, Rineq, k)
%
%   L     : n x n candidate impact matrix, chol(Sigtilde,'lower')*Q for a random
%           rotation Q (see bvar.structural.qr_sign)
%   S     : n x m sign restrictions, one column per shock; +1 and -1 restrict the
%           sign of that response on impact, NaN leaves it free
%   Rineq : m x n x k RANKING restrictions (the paper's term; the legacy code
%           calls them row inequalities), k per shock, each a linear combination
%           of impact responses required to be <= 0 for that shock's column.
%           Pass an m x n matrix when k = 1
%   k     : number of ranking restrictions per shock
%   ok    : true if an admissible assignment exists
%
%   The ranking test here is Rineq*L(:,j) <= 0, so a ROW OF ZEROS is satisfied
%   and imposes nothing. That makes zeros(m,n) the natural way to say "sign
%   restrictions only". bvar.structural.sign_restrict tests the same quantity
%   STRICTLY, < 0, where a zero row instead rejects every candidate and yields an
%   empty identified set with no error. The two are not interchangeable on this
%   point: with no ranking restrictions, pass zeros(m,n) here and an EMPTY Ridx
%   there.
%   L     : on acceptance, the columns reordered so that column i is the shock-i
%           column, with signs flipped where that is what made the restrictions
%           hold, and the remaining n-m columns randomly permuted and re-signed;
%           unchanged when ok is false
%
% WHY IT IS FASTER. bvar.structural.sign_restrict requires column i to satisfy
% shock i and rejects the draw at the first shock that does not, which is the
% standard rejection scheme. But the labelling of the columns of Q is arbitrary:
% a rotation whose third column satisfies the monetary restrictions is just as
% admissible as one where the monetary shock happens to land in the third
% column, and the strict scheme discards the first. This function builds the
% m x n table of which columns admit which shocks, accepts whenever every shock
% has at least one, and draws an assignment uniformly from those available.
%
% It is not merely a heuristic that keeps more draws. Proposition 1 of the paper
% establishes that the accepted R* equals L*Q* for a Q* that is still uniform on
% the orthogonal group, so the target distribution is unchanged; the proof turns
% on the Haar measure being invariant to right multiplication by a permutation
% and a sign matrix. The reported gain at n = 15 with 1000 admissible draws is
% about 3.6 billion candidate rotations and six days for the rejection scheme
% against about 31,000 and sixteen seconds here.
%
% TWO THINGS THE CALLER MUST GET RIGHT, both from the paper. Accept or reject
% the pair (A,Sigma) and Q JOINTLY, as the loop in the replication driver does;
% resampling Q against a fixed posterior draw until it passes targets a
% different distribution (Arias, Rubio-Ramirez, Shin and Waggoner, 2024). And
% take ONE draw per accepted pair, not several: two assignments from the same
% (Sigma,Q) differ only by a permutation and sign flips, so they are dependent.
%
% Both functions are correct and neither replaces the other: sign_restrict is
% the scheme the earlier papers use and the one their replication code
% reproduces, so it stays. See the never-merge list in tests/variant_map.md.
%
% rng consumption on ACCEPTANCE, in order: one unidrnd per shock (m draws,
% choosing among that shock's admissible columns), then randperm(n-m) and
% rand(n-m,1) for the unrestricted columns. A rejected candidate consumes
% nothing, so the stream position depends on how many draws were accepted.
%
% THE CONDITION. Each shock draws its column independently, with no check that a
% column is already taken. What makes that safe is a requirement on the
% restrictions: any two shocks must be separable by their impact responses
% alone, having two common variables on which they agree in sign on one and
% disagree on the other. Every column then admits at most one shock. The paper
% states this as Assumption 1 for sign restrictions and, in the same terms, as
% Assumption 2 once ranking restrictions are included. Restrictions violating it
% are out of scope: the function errors rather than mis-assigns, and the paper's
% second algorithm, which enumerates the admissible set instead, is not
% implemented here.
%
% Body from chan_matthes_yu2026_qe_svarsign/legacy/proposed_15var.m lines 72-113,
% wrapped as a function: m and n come from size(S) and size(L), and the
% acceptance test nnz(sum(abs(satTab),2)) == m is returned as ok rather than
% gating the caller's storage block.
% Equivalence: tests/unit/test_sign_assign.m. Record: tests/variant_map.md.
%
% See:
% Chan, J.C.C., Matthes, C. and Yu, X. (2026). Large Structural VARs with
% Multiple Sign and Ranking Restrictions, Quantitative Economics, 17(3): 709-740.

function [ok,L] = sign_assign(L,S,Rineq,k)
[n,m] = size(S);
if nargin < 4 || isempty(k)
    k = 1;
end

satTab = zeros(m,n); % (i,j) = 1 if the j-th column of L satisfies all restrictions for the i-th shock
                     % (i,j) = -1 if the negative of j-th column of L satisfies all restrictions
for i=1:m  % check sign restrictions & row inequilities
    idx = find(S(:,i)==-1 | S(:,i)==1);
    nidx = length(idx);
    signL = sign(L(idx,:));
    for j=1:n
        if k == 1
            if (sum(signL(:,j) == S(idx,i)) == nidx) && ...
                    (sum(Rineq(i,:,:)*L(:,j) <= 0) == 1)
                satTab(i,j) = 1;
            elseif (sum(signL(:,j) == -S(idx,i)) == nidx) && ...
                    (sum(Rineq(i,:,:)*(-L(:,j)) <= 0) == 1)
                satTab(i,j) = -1;
            end
        elseif k>1
            if (sum(signL(:,j) == S(idx,i)) == nidx) && ...
                    (sum(squeeze(Rineq(i,:,:))'*L(:,j) <= 0) == k)
                satTab(i,j) = 1;
            elseif (sum(signL(:,j) == -S(idx,i)) == nidx) && ...
                    (sum(squeeze(Rineq(i,:,:))'*(-L(:,j)) <= 0) == k)
                satTab(i,j) = -1;
            end
        end
    end
end

ok = (nnz(sum(abs(satTab),2)) == m);  % admissible set is non-empty
if ok
    reorder = zeros(1,n);
    for i=1:m
        idx = find(satTab(i,:));
        draw = idx(unidrnd(length(idx)));
        reorder(i) = draw;
        if satTab(i,draw) == -1
            L(:,draw) = -L(:,draw);
        end
    end
    reorder(m+1:end) = setdiff(1:n,reorder);
    L = L(:,reorder);
        % randomly permute and switch the signs of the last n-m columns
    tmpL2 = L(:,m+1:end);
    tmpL2 = tmpL2(:,randperm(n-m))*diag(2*(rand(n-m,1)>.5) - 1);
    L(:,m+1:end) = tmpL2;
end
end
