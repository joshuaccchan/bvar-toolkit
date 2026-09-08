% This is the main script for replicating the application in section 3 (15-VAR): Figure 1

clear; clc;
rng(123);
% Figure 1(a), RWZ
RWZ_15var; save('results_15var_RWZ.mat','store_response');
horizon = 36;
response_median = squeeze(median(store_response));
response_CI = squeeze(quantile(store_response,[.16,.84]));
ylim_u = [.008, .003, .4, .03, .015, .2];
ylim_l = [-.001, -.001, -.2, -.01, -.005, -.4];
titletext = ["GDP" "GDP Deflator" "3-month Tbill" "Investment" "S&P 500" "Spread"];
figure;
for jj = 1:6
    subplot(3,2,jj);
    hold on
    plotCI((1:horizon-1)',squeeze(response_CI(1,jj,5,2:end)),squeeze(response_CI(2,jj,5,2:end)));
    plot(1:horizon-1,squeeze(response_median(jj,5,2:end)),'--k','LineWidth',1);
    hold off
    line(xlim, [0,0], 'Color', 'k', 'LineWidth', .5); % Draw line for x-axis
    grid on; title(titletext(jj)); xlim([.95, horizon-1]);
    ylim([ylim_l(jj) ylim_u(jj)]); box off;
end
set(gcf,'Position',[100 300 600 400]);
print(gcf, 'Fig1_Panela', '-depsc', '-painters');

clear; clc;
% Figure 1(b), proposed method 
proposed_15var; save('results_15var.mat','store_response');
horizon = 36;
response_median = squeeze(median(store_response));
response_CI = squeeze(quantile(store_response,[.16,.84]));
ylim_u = [.008, .003, .4, .03, .015, .2];
ylim_l = [-.001, -.001, -.2, -.01, -.005, -.4];
titletext = ["GDP" "GDP Deflator" "3-month Tbill" "Investment" "S&P 500" "Spread"];
figure;
for jj = 1:6
    subplot(3,2,jj);
    hold on
    plotCI((1:horizon-1)',squeeze(response_CI(1,jj,5,2:end)),squeeze(response_CI(2,jj,5,2:end)));
    plot(1:horizon-1,squeeze(response_median(jj,5,2:end)),'--k','LineWidth',1);
    hold off
    line(xlim, [0,0], 'Color', 'k', 'LineWidth', .5); % Draw line for x-axis
    grid on; title(titletext(jj)); xlim([.95, horizon-1]);
    ylim([ylim_l(jj) ylim_u(jj)]); box off;
end
set(gcf,'Position',[100 300 600 400]);
print(gcf, 'Fig1_Panelb', '-depsc', '-painters');


