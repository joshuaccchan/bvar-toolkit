% bvar.forecast.tables - forecast-evaluation accumulation and RMSFE / average
% log predictive likelihood (ALPL) table construction, replicated VERBATIM from
% the recursive-forecasting driver tails. Extracted 2026-09-01 (step 6,
% forecast engine). Three named actions:
%
%   row = bvar.forecast.tables('accum_row', tmpyhat, obs)
%     One accumulation row for one vintage and one horizon:
%       tmpmax = max(tmpyhat(:,n+1:end));
%       row = [obs mean(tmpyhat(:,1:n)) log(mean(exp(tmpyhat(:,n+1:end)-tmpmax)))+tmpmax];
%     with n = (size(tmpyhat,2)-1)/2. tmpyhat is the nsim x (2n+1) per-draw
%     matrix from bvar.forecast.iterate ([point forecasts, per-variable log
%     prelikes, joint log prelike]); obs is the 1 x n observed outturn row.
%     Canonicalizes chan2020_springer_largebvar/legacy/main_forecasting.m
%     lines 147-154 (yhat0 and yhat1 rows) and chan2021_ijf_mahp/legacy/
%     main_forecasting.m lines 98-105 (yhat1 and yhat4 rows) - the same
%     formula in both packages. The storage GUARDS (springer: yhat1 only when
%     t<=T-1; MAHP: yhat4 only when t<=T-4) are loop control and stay with the
%     caller. NOTE on older MATLABs tmpyhat can arrive complex-typed with zero
%     imaginary part (see bvar.forecast.iterate header; on R2025b diag() demotes
%     it back to real) and max() then compares by magnitude - the formula is
%     reproduced verbatim so this behavior is preserved exactly either way.
%     NOT covered (deferred to the kronecker family pass): the subsetted
%     accumulation of chan2020_jbes_kronecker/legacy/main_forecasting.m
%     (columns [n+var_small, 2*n+1] only, 3*4+1-wide rows).
%
%   S = bvar.forecast.tables('springer', yhat0, yhat1, var_core)
%     RMSFE / ALPL tables of chan2020_springer_largebvar/legacy/
%     main_forecasting.m lines 160-172, evaluation-period row trims verbatim
%     (yhat0(5:end,:) = 1975Q1 on; yhat1(4:end,:)). var_core nonempty (the
%     legacy [1 7 8 12]' column) gives the model~=1 form: RMSFE over the
%     var_core columns and ALPL over columns 2*n+var_core (per-variable only -
%     the joint column is NOT tabled). var_core = [] gives the model==1
%     (BVAR-small) form: all n variables, ALPL over 2*n+1:end (per-variable
%     AND joint). S has fields RMSFE_0, RMSFE_1, RMSFE = [RMSFE_0 RMSFE_1],
%     aveprelike_0, aveprelike_1, aveprelike.
%
%   S = bvar.forecast.tables('mahp', yhat1, yhat4)
%     RMSFE / ALPL tables of chan2021_ijf_mahp/legacy/main_forecasting.m
%     lines 111-116, row trims verbatim (yhat1(4:end,:), yhat4(1:end,:)),
%     all n variables, ALPL over 2*n+1:end (per-variable AND joint). S has
%     fields RMSFE_1, RMSFE_4, RMSFE, aveprelike_1, aveprelike_4, aveprelike.
%
% The legacy fprintf display blocks (headline-variable pretty-printing) are
% formatting only and are not reproduced; callers print from S.
%
% See:
% Chan, J.C.C. (2020). Large Bayesian Vector Autoregressions. In: P. Fuleky (Ed.),
% Macroeconomic Forecasting in the Era of Big Data, 95-125, Springer, Cham.
% Chan, J.C.C. (2021). Minnesota-Type Adaptive Hierarchical Priors for
% Large Bayesian VARs, International Journal of Forecasting, 37(3), 1212-1226.

function out = tables(action, varargin)
switch action
    case 'accum_row'
        tmpyhat = varargin{1}; obs = varargin{2};
        n = (size(tmpyhat,2)-1)/2;
        tmpmax = max(tmpyhat(:,n+1:end));
        out = [obs mean(tmpyhat(:,1:n)) ...
            log(mean(exp(tmpyhat(:,n+1:end)-tmpmax)))+tmpmax];

    case 'springer'
        yhat0 = varargin{1}; yhat1 = varargin{2}; var_core = varargin{3};
        n = (size(yhat0,2)-1)/3;
        if isempty(var_core)    % model 1 (BVAR-small): all n variables
            RMSFE_0 = sqrt(mean((yhat0(5:end,1:n) - yhat0(5:end,n+1:2*n)).^2))'; % evaluation period: 1975Q1-2015Q4
            RMSFE_1 = sqrt(mean((yhat1(4:end,1:n) - yhat1(4:end,n+1:2*n)).^2))';
            aveprelike_0 = mean(yhat0(5:end,2*n+1:end))';
            aveprelike_1 = mean(yhat1(4:end,2*n+1:end))';
        else                    % models 2-8: headline var_core columns
            RMSFE_0 = sqrt(mean((yhat0(5:end,var_core) - yhat0(5:end,n+var_core)).^2))'; % evaluation period: 1975Q1-2015Q4
            RMSFE_1 = sqrt(mean((yhat1(4:end,var_core) - yhat1(4:end,n+var_core)).^2))';
            aveprelike_0 = mean(yhat0(5:end,2*n+var_core))';
            aveprelike_1 = mean(yhat1(4:end,2*n+var_core))';
        end
        out = struct('RMSFE_0',RMSFE_0,'RMSFE_1',RMSFE_1, ...
            'RMSFE',[RMSFE_0 RMSFE_1], ...
            'aveprelike_0',aveprelike_0,'aveprelike_1',aveprelike_1, ...
            'aveprelike',[aveprelike_0 aveprelike_1]);

    case 'mahp'
        yhat1 = varargin{1}; yhat4 = varargin{2};
        n = (size(yhat1,2)-1)/3;
        RMSFE_1 = sqrt(mean((yhat1(4:end,1:n) - yhat1(4:end,n+1:2*n)).^2))';
        RMSFE_4 = sqrt(mean((yhat4(1:end,1:n) - yhat4(1:end,n+1:2*n)).^2))';
        aveprelike_1 = mean(yhat1(4:end,2*n+1:end))';
        aveprelike_4 = mean(yhat4(1:end,2*n+1:end))';
        out = struct('RMSFE_1',RMSFE_1,'RMSFE_4',RMSFE_4, ...
            'RMSFE',[RMSFE_1 RMSFE_4], ...
            'aveprelike_1',aveprelike_1,'aveprelike_4',aveprelike_4, ...
            'aveprelike',[aveprelike_1 aveprelike_4]);

    otherwise
        error('bvar:forecast:tables:unknownAction', ...
            'unknown action ''%s''; use accum_row, springer or mahp', action);
end
end
