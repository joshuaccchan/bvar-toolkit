% This is the main script for replicating the application in Uhlig (2005):
% Table 2 and Figures 11-13 (Appendix B）
% NOTICE: In compliance with licensing agreements, the commodity prices 
% presented have been perturbed with noise; the empirical results are virtually unaffected.

clear; clc;
rng(123);
% Figure 11, RWZ
algom = 1;     % 1: standard accept-reject algorithm; 2: proposed algorithm; 3: Read(2022) algorithm
Application_Uhlig2005;
time_cost1 = etime(clock,start_time)/60;
horizon = 61;
response_median = squeeze(median(store_response));
response_CI = squeeze(quantile(store_response,[.16,.84]));
titletext = ["output" "inflation" "commodity prices" ...
    "nonborrowed reserve" "total reserve" "FFR"];
pos = [1,3,5,4,2,6];
figure;
for jj = 1:6
    subplot(3,2,pos(jj));
    hold on
    plotCI((1:horizon-1)',squeeze(response_CI(1,jj,2:end)),squeeze(response_CI(2,jj,2:end)));
    plot(1:horizon-1,squeeze(response_median(jj,2:end)),'--k','LineWidth',1);
    hold off
    line(xlim, [0,0], 'Color', 'k', 'LineWidth', .5); % Draw line for x-axis
    grid on; title(titletext(jj)); xlim([.95, horizon-1]);
    box off;
end
set(gcf,'Position',[100 300 600 400]);
print(gcf, 'Fig11', '-depsc', '-painters');

% Figure 12, proposed algorithm
algom = 2;     % 1: standard accept-reject algorithm; 2: proposed algorithm; 3: Read(2022) algorithm
Application_Uhlig2005;
time_cost2 = etime(clock,start_time)/60;
horizon = 61;
response_median = squeeze(median(store_response));
response_CI = squeeze(quantile(store_response,[.16,.84]));
titletext = ["output" "inflation" "commodity prices" ...
    "nonborrowed reserve" "total reserve" "FFR"];
pos = [1,3,5,4,2,6];
figure;
for jj = 1:6
    subplot(3,2,pos(jj));
    hold on
    plotCI((1:horizon-1)',squeeze(response_CI(1,jj,2:end)),squeeze(response_CI(2,jj,2:end)));
    plot(1:horizon-1,squeeze(response_median(jj,2:end)),'--k','LineWidth',1);
    hold off
    line(xlim, [0,0], 'Color', 'k', 'LineWidth', .5); % Draw line for x-axis
    grid on; title(titletext(jj)); xlim([.95, horizon-1]);
    box off;
end
set(gcf,'Position',[100 300 600 400]);
print(gcf, 'Fig12', '-depsc', '-painters');

% Figure 13, Read(2022) algorithm
algom = 3;     % 1: standard accept-reject algorithm; 2: proposed algorithm; 3: Read(2022) algorithm
Application_Uhlig2005;
time_cost3 = etime(clock,start_time)/60;
horizon = 61;
response_median = squeeze(median(store_response));
response_CI = squeeze(quantile(store_response,[.16,.84]));
titletext = ["output" "inflation" "commodity prices" ...
    "nonborrowed reserve" "total reserve" "FFR"];
pos = [1,3,5,4,2,6];
figure;
for jj = 1:6
    subplot(3,2,pos(jj));
    hold on
    plotCI((1:horizon-1)',squeeze(response_CI(1,jj,2:end)),squeeze(response_CI(2,jj,2:end)));
    plot(1:horizon-1,squeeze(response_median(jj,2:end)),'--k','LineWidth',1);
    hold off
    line(xlim, [0,0], 'Color', 'k', 'LineWidth', .5); % Draw line for x-axis
    grid on; title(titletext(jj)); xlim([.95, horizon-1]);
    box off;
end
set(gcf,'Position',[100 300 600 400]);
print(gcf, 'Fig13', '-depsc', '-painters');

% Table 2, time cost comparison
Table2 = [time_cost2, time_cost1, time_cost3] %proposed method, RWZ, Read
csvwrite('Table2.csv',Table2)