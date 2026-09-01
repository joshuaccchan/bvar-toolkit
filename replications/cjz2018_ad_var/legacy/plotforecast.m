% Chan, J.C.C., Jacobi, L. and Zhu, D. (2018). How Sensitive Are VAR 
% Forecasts to Prior Hyperparameters? An Automated Sensitivity Analysis,
% CAMA Working Paper 25/2018.

%% forecat plot
t=T-19:T+forecast_period;
if is_fullsample
    tid = linspace(1955,2017.75+forecast_period/4,t(end))';
else
    tid = linspace(2000,2017.75+forecast_period/4,t(end))';
end
figure; 
subplot(3,1,1);hold on;
plot(tid(T-19:T),Y(T-19:end,1),'k','linewidth',1,'DisplayName','data');
plot(tid(T:end),[Y(end,1);means(:,1)],'b-','linewidth',1,'DisplayName','point');
plot(tid(T:end),[Y(end,1);Q_upper(:,1)],'r--','linewidth',1,'DisplayName','84%');
plot(tid(T:end),[Y(end,1);Q_lower(:,1)],'r:','linewidth',1,'DisplayName','16%');
ylim([3 8]); box off; title('Unemployment'); 
% legend('show');
hold off;

subplot(3,1,2);hold on;
plot(tid(T-19:T),Y(T-19:end,2),'k','linewidth',1,'DisplayName','data');
plot(tid(T:end),[Y(end,2);means(:,2)],'b-','linewidth',1,'DisplayName','point');
plot(tid(T:end),[Y(end,2);Q_upper(:,2)],'r--','linewidth',1,'DisplayName','84%');
plot(tid(T:end),[Y(end,2);Q_lower(:,2)],'r:','linewidth',1,'DisplayName','16%');
ylim([-3 9]); box off; title('Interest Rate'); 
% legend('show');
hold off;

subplot(3,1,3);hold on;
plot(tid(T-19:T),Y(T-19:end,3),'k','linewidth',1,'DisplayName','data');
plot(tid(T:end),[Y(end,3);means(:,3)],'b-','linewidth',1,'DisplayName','point');
plot(tid(T:end),[Y(end,3);Q_upper(:,3)],'r--','linewidth',1,'DisplayName','84%');
plot(tid(T:end),[Y(end,3);Q_lower(:,3)],'r:','linewidth',1,'DisplayName','16%');
ylim([-4 10]); box off; title('GDP Growth'); 
legend('show');
hold off;


%% Sensitivity plot
% subplot(3,3,1);
% plot(1:forecast_period, squeeze(Ym.d(:,:,1))');
% subplot(3,3,2);
% plot(1:forecast_period, squeeze(Ym.d(:,:,2))');
% subplot(3,3,3);
% plot(1:forecast_period, squeeze(Ym.d(:,:,3))');
% 
% subplot(3,3,4);
% plot(1:forecast_period, squeeze(d_Qupper(:,:,1))');
% subplot(3,3,5);
% plot(1:forecast_period, squeeze(d_Qupper(:,:,2))');
% subplot(3,3,6);
% plot(1:forecast_period, squeeze(d_Qupper(:,:,3))');
% 
% subplot(3,3,7);
% plot(1:forecast_period, squeeze(d_QLower(:,:,1))');
% subplot(3,3,8);
% plot(1:forecast_period, squeeze(d_QLower(:,:,2))');
% subplot(3,3,9);
% plot(1:forecast_period, squeeze(d_QLower(:,:,3))');

figure;
subplot(3,3,1);
hold on;
plot(1:forecast_period, Ym.d(:,1,1)','b-','linewidth',1);
plot(1:forecast_period, Ym.d(:,2,1)','r--','linewidth',1);
plot(1:forecast_period, Ym.d(:,3,1)','k-.','linewidth',1);
hold off;
box off; title('\kappa_1'); ylabel('point');

subplot(3,3,2);
hold on;
plot(1:forecast_period, Ym.d(:,1,2)','b-','linewidth',1);
plot(1:forecast_period, Ym.d(:,2,2)','r--','linewidth',1);
plot(1:forecast_period, Ym.d(:,3,2)','k-.','linewidth',1);
hold off;
box off; title('\kappa_2');

subplot(3,3,3);
hold on;
plot(1:forecast_period, Ym.d(:,1,3)','b-','linewidth',1,'DisplayName','GDP');
plot(1:forecast_period, Ym.d(:,2,3)','r--','linewidth',1,'DisplayName','i');
plot(1:forecast_period, Ym.d(:,3,3)','k-.','linewidth',1,'DisplayName','u');
hold off;
box off; title('\kappa_3');

subplot(3,3,4);
hold on;
plot(1:forecast_period, d_Qupper(1,:,1)','b-','linewidth',1);
plot(1:forecast_period, d_Qupper(2,:,1)','r--','linewidth',1);
plot(1:forecast_period, d_Qupper(3,:,1)','k-.','linewidth',1);
hold off;
box off; ylabel('84%');

subplot(3,3,5);
hold on;
plot(1:forecast_period, d_Qupper(1,:,2)','b-','linewidth',1);
plot(1:forecast_period, d_Qupper(2,:,2)','r--','linewidth',1);
plot(1:forecast_period, d_Qupper(3,:,2)','k-.','linewidth',1);
hold off;
box off; 

subplot(3,3,6);
hold on;
plot(1:forecast_period, d_Qupper(1,:,3)','b-','linewidth',1,'DisplayName','GDP');
plot(1:forecast_period, d_Qupper(2,:,3)','r--','linewidth',1,'DisplayName','i');
plot(1:forecast_period, d_Qupper(3,:,3)','k-.','linewidth',1,'DisplayName','u');
hold off;
box off;

subplot(3,3,7);
hold on;
plot(1:forecast_period, d_QLower(1,:,1)','b-','linewidth',1);
plot(1:forecast_period, d_QLower(2,:,1)','r--','linewidth',1);
plot(1:forecast_period, d_QLower(3,:,1)','k-.','linewidth',1);
hold off;
box off; ylabel('16%');

subplot(3,3,8);
hold on;
plot(1:forecast_period, d_QLower(1,:,2)','b-','linewidth',1);
plot(1:forecast_period, d_QLower(2,:,2)','r--','linewidth',1);
plot(1:forecast_period, d_QLower(3,:,2)','k-.','linewidth',1);
hold off;
box off; 

subplot(3,3,9);
hold on;
plot(1:forecast_period, d_QLower(1,:,3)','b-','linewidth',1,'DisplayName','GDP');
plot(1:forecast_period, d_QLower(2,:,3)','r--','linewidth',1,'DisplayName','i');
plot(1:forecast_period, d_QLower(3,:,3)','k-.','linewidth',1,'DisplayName','u');
hold off;
box off; legend('show');

