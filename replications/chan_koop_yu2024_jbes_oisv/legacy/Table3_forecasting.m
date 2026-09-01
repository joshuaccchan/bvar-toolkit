n = 20;
addpath('./results_mat')
M1 = load('forecastingCS1-cluster.mat'); M1_yhat1 = M1.yhatmedian(:,:,1); M1_yhat6 = M1.yhatmedian(1:end-5,:,6); M1_yhat12 = M1.yhatmedian(1:end-11,:,12);
M2 = load('forecastingCS2-cluster.mat'); M2_yhat1 = M2.yhatmedian(:,:,1); M2_yhat6 = M2.yhatmedian(1:end-5,:,6); M2_yhat12 = M2.yhatmedian(1:end-11,:,12);
M3 = load('forecastingOI1-cluster.mat'); M3_yhat1 = M3.yhatmedian(:,:,1); M3_yhat6 = M3.yhatmedian(1:end-5,:,6); M3_yhat12 = M3.yhatmedian(1:end-11,:,12);
M4 = load('forecastingOI2-cluster.mat'); M4_yhat1 = M4.yhatmedian(:,:,1); M4_yhat6 = M4.yhatmedian(1:end-5,:,6); M4_yhat12 = M4.yhatmedian(1:end-11,:,12);

T = length(M1_yhat12(:,1));
Table_re = []; Table_DM = []; DMbenchmark = 1;
for i = 1:4
    
    R1 = zeros(T,4); DM_R1 = zeros(4,1);
    R1(:,1) = M1_yhat1(12:end,i) - M1_yhat1(12:end,n+i);
    R1(:,2) = M2_yhat1(12:end,n+1-i) - M2_yhat1(12:end,2*n+1-i); % for reversed model, same for m4
    R1(:,3) = M3_yhat1(12:end,i) - M3_yhat1(12:end,n+i);
    R1(:,4) = M4_yhat1(12:end,n+1-i) - M4_yhat1(12:end,2*n+1-i);
    for ii = 1:4
        if ii ~= DMbenchmark
            tempre = EvalFore([R1(:,DMbenchmark), R1(:,ii)],2); DM_R1(ii) = [tempre.DM(2)];
        else
            DM_R1(ii) = NaN;
        end
    end
    
    R6 = zeros(T,4); DM_R6 = zeros(4,1);
    R6(:,1) = M1_yhat6(7:end,i) - M1_yhat6(7:end,n+i);
    R6(:,2) = M2_yhat6(7:end,n+1-i) - M2_yhat6(7:end,2*n+1-i);
    R6(:,3) = M3_yhat6(7:end,i) - M3_yhat6(7:end,n+i);
    R6(:,4) = M4_yhat6(7:end,n+1-i) - M4_yhat6(7:end,2*n+1-i);
    for ii = 1:4
        if ii ~= DMbenchmark
            tempre = EvalFore([R6(:,DMbenchmark), R6(:,ii)],2); DM_R6(ii) = [tempre.DM(2)];
        else
            DM_R6(ii) = NaN;
        end
    end
    
    R12 = zeros(T,4); DM_R12 = zeros(4,1);
    R12(:,1) = M1_yhat12(1:end,i) - M1_yhat12(1:end,n+i);
    R12(:,2) = M2_yhat12(1:end,n+1-i) - M2_yhat12(1:end,2*n+1-i);
    R12(:,3) = M3_yhat12(1:end,i) - M3_yhat12(1:end,n+i);
    R12(:,4) = M4_yhat12(1:end,n+1-i) - M4_yhat12(1:end,2*n+1-i);
    for ii = 1:4
        if ii ~= DMbenchmark
            tempre = EvalFore([R12(:,DMbenchmark), R12(:,ii)],2); DM_R12(ii) = [tempre.DM(2)];
        else
            DM_R12(ii) = NaN;
        end
    end
    
    L1 = zeros(T,4); DM_L1 = zeros(4,1);
    L1(:,1) = M1_yhat1(12:end,2*n+i);
    L1(:,2) = M2_yhat1(12:end,3*n+1-i);
    L1(:,3) = M3_yhat1(12:end,2*n+i);
    L1(:,4) = M4_yhat1(12:end,3*n+1-i);
    for ii = 1:4
        if ii ~= DMbenchmark
            tempre = EvalFore([L1(:,DMbenchmark), L1(:,ii)],2); DM_L1(ii) = [tempre.DM(2)];
        else
            DM_L1(ii) = NaN;
        end
    end
    
    L6 = zeros(T,4); DM_L6 = zeros(4,1);
    L6(:,1) = M1_yhat6(7:end,2*n+i);
    L6(:,2) = M2_yhat6(7:end,3*n+1-i);
    L6(:,3) = M3_yhat6(7:end,2*n+i);
    L6(:,4) = M4_yhat6(7:end,3*n+1-i);
    for ii = 1:4
        if ii ~= DMbenchmark
            tempre = EvalFore([L6(:,DMbenchmark), L6(:,ii)],2); DM_L6(ii) = [tempre.DM(2)];
        else
            DM_L6(ii) = NaN;
        end
    end
    
    L12 = zeros(T,4); DM_L12 = zeros(4,1);
    L12(:,1) = M1_yhat12(1:end,2*n+i);
    L12(:,2) = M2_yhat12(1:end,3*n+1-i);
    L12(:,3) = M3_yhat12(1:end,2*n+i);
    L12(:,4) = M4_yhat12(1:end,3*n+1-i);
    for ii = 1:4
        if ii ~= DMbenchmark
            tempre = EvalFore([L12(:,DMbenchmark), L12(:,ii)],2); DM_L12(ii) = [tempre.DM(2)];
        else
            DM_L12(ii) = NaN;
        end
    end
    
    Table_re = [Table_re; [sqrt(mean(R1.^2))', sqrt(mean(R6.^2))', sqrt(mean(R12.^2))', mean(L1)', mean(L6)', mean(L12)']];
    Table_DM = [Table_DM; [DM_R1, DM_R6, DM_R12, DM_L1, DM_L6, DM_L12]];
end
Table_DM(Table_DM<=0.01)=3; Table_DM(Table_DM<=0.05)=2; Table_DM(Table_DM<=0.10)=1; Table_DM(Table_DM<1)=0; 

Table_re_withDM = [];
for i = 1:6
    tempre = [];
    if i <=3 option = 1; else option = 2; end % option-1: min is optimum; option-2: max is optimum;
    for j = 1:4
        tempre = [tempre; print_TexRe(Table_re(4*(j-1)+1:4*j,i),Table_DM(4*(j-1)+1:4*j,i),option)];
    end
    Table_re_withDM = [Table_re_withDM, tempre];
end

Table_re_withDM = cell2table(Table_re_withDM)
writetable(Table_re_withDM, strcat('Table3_LaTex_cluster.xls'),'WriteVariableNames',0);