% bvar.util.gam_mode - the modal configuration of a stored binary indicator
% matrix: the row value that occurs most often across draws. Used on the hybrid
% TVP-VAR's store_gam, whose nsim x 2n rows record which equations had
% time-varying VAR coefficients and which had time-varying impact elements, to
% report the single most-visited model rather than a marginal mean.
%
%   gam_mode = bvar.util.gam_mode(store_gam)
%
%   store_gam : nsim x nmodel matrix of 0/1 draws
%   gam_mode  : 1 x nmodel row, the most frequent configuration
%
% Each row is keyed by its binary expansion gam*2.^(0:nmodel-1)', counts are
% accumulated by growing a table, and ties are broken by whichever configuration
% sortrows leaves on top. The loop over draws is the legacy implementation and is
% O(nsim x distinct configurations); it is quick at the sizes used here but slows
% down for a long chain over many equations.
%
% Body from chan2023_jbes_hybtvp/legacy/utility/get_gammode.m, renamed.
% Equivalence: tests/unit/test_hybtvp_equivalence.m. Record: tests/variant_map.md.
%
% See:
% Chan, J.C.C. (2023). Large Hybrid Time-Varying Parameter VARs, Journal of
% Business and Economic Statistics, 41(3): 890-905.

function gam_mode = gam_mode(store_gam)
[nsim,nmodel] = size(store_gam);
model_idx = store_gam(1,:)*(2.^(0:nmodel-1))' + 1;
model_count = [store_gam(1,:), model_idx, 1]; % [gam, model id, count]
for ii=2:nsim
    gam = store_gam(ii,:);
    model_idx = gam*(2.^(0:nmodel-1))' + 1;
    id = find(model_count(:,nmodel+1) == model_idx);
    if isempty(id)
        model_count = [model_count; [gam,model_idx,1]]; %#ok<AGROW> % legacy growth pattern, kept verbatim
    else
        model_count(id,end) = model_count(id,end)+1;
    end
end
model_count = sortrows(model_count,-(nmodel+2));
gam_mode = model_count(1,1:nmodel);
end
