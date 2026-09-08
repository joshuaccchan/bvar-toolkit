function test_sign_assign
% bvar.structural.sign_assign must reproduce the inline block of
% proposed_15var.m draw for draw, on the paper's own 15-variable restrictions,
% and must leave the rng in the same state. Also checks the two structural
% properties the function relies on: that an accepted L has each shock's
% restrictions satisfied by ITS OWN column after the reordering, and that
% sign_assign accepts strictly more often than bvar.structural.sign_restrict,
% which is the reason it exists.
root = getappdata(0, 'bvar_repo_root');
leg = fullfile(root, 'replications', 'chan_matthes_yu2026_qe_svarsign', 'legacy');
addpath(fullfile(leg, 'utility'));
c = onCleanup(@() rmpath(fullfile(leg, 'utility'))); %#ok<NASGU>

data = readmatrix(fullfile(leg,'data','database_2019Q4.csv'),'Range','B3:P142');
p = 5; idx_ns = [1,2,4,5,10,11,12,13,15];
Y0 = data(1:8,1:15); Y = data(9:end,1:15);
[~,n] = size(Y);
sig2 = bvar.priors.resid_var_ar4(Y0,Y);
prior_redu = bvar.priors.acp_redu(n,p,[1,1,1,100],sig2,idx_ns);

S = [[1,-1,NaN,NaN,1,NaN,NaN,NaN,NaN,-1,-1,NaN,1,NaN,1]', ...
     [1,1,1,NaN,NaN,NaN,NaN,NaN,NaN,1,1,NaN,1,1,NaN]', ...
     [1,1,-1,NaN,NaN,NaN,NaN,NaN,NaN,1,1,NaN,1,-1,NaN]', ...
     [1,1,1,NaN,-1,NaN,NaN,NaN,NaN,1,1,NaN,1,1,-1]', ...
     [1,1,1,NaN,1,NaN,NaN,NaN,NaN,1,1,NaN,1,1,1]'];
m = size(S,2);
Rineq = [zeros(1,n); -1,0,0,1,zeros(1,11); zeros(1,n); ...
         1,0,0,-1,zeros(1,11); 1,0,0,-1,zeros(1,11)];

nbatch = 400;
rng(20260908,'twister');
[sa,sb,ss] = bvar.samplers.acp_theta_sig(Y0,Y,p,prior_redu,nbatch);
[~,store_Sigtilde] = bvar.structural.reduced_form(sa,sb,ss);

nacc = 0; nacc_strict = 0;
for isim = 1:nbatch
    Sigtilde = squeeze(store_Sigtilde(isim,:,:));
    [Q,~] = bvar.structural.qr_sign(randn(n,n));
    L = chol(Sigtilde,'lower')*Q;

        % draw-for-draw against the legacy block, from the same rng state
    s0 = rng;
    [ok, Lout] = bvar.structural.sign_assign(L, S, Rineq, 1);
    s1 = rng;
    rng(s0);
    [okref, Lref] = legacy_block(L, S, Rineq, 1);
    assert(isequal(ok, okref), 'draw %d: acceptance differs', isim);
    assert(isequal(Lout, Lref), 'draw %d: L differs', isim);
    assert(isequal(s1.State, subsref(rng, substruct('.','State'))), ...
        'draw %d: rng consumption differs', isim);

    if ok
        nacc = nacc + 1;
            % after reordering, column i must satisfy shock i's sign pattern
        for i = 1:m
            idx = find(S(:,i)==-1 | S(:,i)==1);
            assert(all(sign(Lout(idx,i)) == S(idx,i)), ...
                'draw %d: column %d does not satisfy its shock after reordering', isim, i);
        end
    end

        % the strict scheme, for the acceptance comparison. Its row inequalities
        % are indexed by column rather than by shock, so the equivalent input is
        % the three nonzero rows applied to the columns they restrict.
    [okS,~] = bvar.structural.sign_restrict(L, S, Rineq([2 4 5],:), [2 4 5]);
    nacc_strict = nacc_strict + okS;
end

assert(nacc > 0, 'no candidate was accepted - the comparison proves nothing');
assert(nacc > nacc_strict, ...
    'sign_assign (%d) did not accept more often than sign_restrict (%d)', nacc, nacc_strict);
end

% -------------------------------------------------------------------------
function [ok,L] = legacy_block(L,S,Rineq,k)
% proposed_15var.m lines 72-113, inline, as the reference
[n,m] = size(S); %#ok<ASGLU>
m = size(S,2); n = size(L,1);
satTab = zeros(m,n);
for i=1:m
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
        end
    end
end
ok = (nnz(sum(abs(satTab),2)) == m);
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
    tmpL2 = L(:,m+1:end);
    tmpL2 = tmpL2(:,randperm(n-m))*diag(2*(rand(n-m,1)>.5) - 1);
    L(:,m+1:end) = tmpL2;
end
end
