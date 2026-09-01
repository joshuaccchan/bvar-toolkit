% Support function for estimating the hybrid TVP-VAR in Chan (2022)
%
% See:
% Chan, J.C.C. (2022). Large Hybrid Time-Varying Parameter VARs, Journal of
% Business and Economic Statistics, forthcoming

function gam_mode = get_gammode(store_gam)
[nsim,nmodel] = size(store_gam);
model_idx = store_gam(1,:)*(2.^(0:nmodel-1))' + 1;
model_count = [store_gam(1,:), model_idx, 1]; % [gam, model id, count]
for ii=2:nsim
    gam = store_gam(ii,:);
    model_idx = gam*(2.^(0:nmodel-1))' + 1;
    id = find(model_count(:,nmodel+1) == model_idx);
    if isempty(id)
        model_count = [model_count; [gam,model_idx,1]];
    else
        model_count(id,end) = model_count(id,end)+1;
    end
end
model_count = sortrows(model_count,-(nmodel+2));
gam_mode = model_count(1,1:nmodel);
end