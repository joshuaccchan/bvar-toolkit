clear all; clc;
addpath('./results_mat')
n = 20;
load('fullsampleOI1.mat'); model=1; rev_option=0;
[var_4_OI1, cov_4_OI1] = get_varcov_first2last2(Sig_mean,rev_option);
load('fullsampleOI2.mat'); model=1; rev_option=1;
[var_4_OI2, cov_4_OI2] = get_varcov_first2last2(Sig_mean,rev_option);
load('fullsampleCS1.mat'); model=2; rev_option=0;
[var_4_CS1, cov_4_CS1] = get_varcov_first2last2(Sig_mean,rev_option);
load('fullsampleCS2.mat'); model=2; rev_option=1;
[var_4_CS2, cov_4_CS2] = get_varcov_first2last2(Sig_mean,rev_option);

startDate = datetime('March 1, 1961');
d1 = startDate + calmonths(0:59*12-3);

%%%% variances:
figure;
subplot(4,1,1); hold on; box off;
plot(d1,var_4_CS1(:,1),'-','color',[0,0,1],'LineWidth',1)
plot(d1,var_4_CS2(:,1),'--','color',[1,0,0],'LineWidth',1)
plot(d1,var_4_OI1(:,1),'-.','color',[0,0,0],'LineWidth',1)
plot(d1,var_4_OI2(:,1),':','color',[0,0,0],'LineWidth',1)
xlabel('variable 1','FontSize',10); %xlim([1 12]);

subplot(4,1,2); hold on; box off;
plot(d1,var_4_CS1(:,2),'-','color',[0,0,1],'LineWidth',1)
plot(d1,var_4_CS2(:,2),'--','color',[1,0,0],'LineWidth',1)
plot(d1,var_4_OI1(:,2),'-.','color',[0,0,0],'LineWidth',1)
plot(d1,var_4_OI2(:,2),':','color',[0,0,0],'LineWidth',1)
xlabel('variable 2','FontSize',10); %xlim([1 12]);

subplot(4,1,3); hold on; box off;
plot(d1,var_4_CS1(:,3),'-','color',[0,0,1],'LineWidth',1)
plot(d1,var_4_CS2(:,3),'--','color',[1,0,0],'LineWidth',1)
plot(d1,var_4_OI1(:,3),'-.','color',[0,0,0],'LineWidth',1)
plot(d1,var_4_OI2(:,3),':','color',[0,0,0],'LineWidth',1)
h=legend('CS-VAR-SV-1','CS-VAR-SV-2','OI-VAR-SV-1','OI-VAR-SV-2');
h.NumColumns=1;
xlabel('variable 19','FontSize',10); %xlim([1 12]);

subplot(4,1,4); hold on; box off;
plot(d1,var_4_CS1(:,4),'-','color',[0,0,1],'LineWidth',1)
plot(d1,var_4_CS2(:,4),'--','color',[1,0,0],'LineWidth',1)
plot(d1,var_4_OI1(:,4),'-.','color',[0,0,0],'LineWidth',1)
plot(d1,var_4_OI2(:,4),':','color',[0,0,0],'LineWidth',1)
xlabel('variable 20','FontSize',10); %xlim([1 12]);

set(gcf, 'Units', 'centimeters', 'Position', [2, 2, 18, 20]);
saveas(gcf,'Fig-var-foursubfigs-first2last2-1960','epsc');
%savefig('Fig-var-foursubfigs-first2last2-1960.fig')

%%%% covariances:
figure;
subplot(4,1,1); hold on; box off;
plot(d1,cov_4_CS1(:,1),'-','color',[0,0,1],'LineWidth',1)
plot(d1,cov_4_CS2(:,1),'--','color',[1,0,0],'LineWidth',1)
plot(d1,cov_4_OI1(:,1),'-.','color',[0,0,0],'LineWidth',1)
plot(d1,cov_4_OI2(:,1),':','color',[0,0,0],'LineWidth',1)
xlabel('covariance betweent variables 1 and 2','FontSize',10); %xlim([1 12]);

subplot(4,1,2); hold on; box off;
plot(d1,cov_4_CS1(:,2),'-','color',[0,0,1],'LineWidth',1)
plot(d1,cov_4_CS2(:,2),'--','color',[1,0,0],'LineWidth',1)
plot(d1,cov_4_OI1(:,2),'-.','color',[0,0,0],'LineWidth',1)
plot(d1,cov_4_OI2(:,2),':','color',[0,0,0],'LineWidth',1)
xlabel('covariance betweent variables 1 and 19','FontSize',10); %xlim([1 12]);

subplot(4,1,3); hold on; box off;
plot(d1,cov_4_CS1(:,3),'-','color',[0,0,1],'LineWidth',1)
plot(d1,cov_4_CS2(:,3),'--','color',[1,0,0],'LineWidth',1)
plot(d1,cov_4_OI1(:,3),'-.','color',[0,0,0],'LineWidth',1)
plot(d1,cov_4_OI2(:,3),':','color',[0,0,0],'LineWidth',1)
h=legend('CS-VAR-SV-1','CS-VAR-SV-2','OI-VAR-SV-1','OI-VAR-SV-2');
h.NumColumns=1;
xlabel('covariance betweent variables 2 and 20','FontSize',10); %xlim([1 12]);

subplot(4,1,4); hold on; box off;
plot(d1,cov_4_CS1(:,4),'-','color',[0,0,1],'LineWidth',1)
plot(d1,cov_4_CS2(:,4),'--','color',[1,0,0],'LineWidth',1)
plot(d1,cov_4_OI1(:,4),'-.','color',[0,0,0],'LineWidth',1)
plot(d1,cov_4_OI2(:,4),':','color',[0,0,0],'LineWidth',1)
xlabel('covariance betweent variables 19 and 20','FontSize',10); %xlim([1 12]);

set(gcf, 'Units', 'centimeters', 'Position', [2, 2, 18, 20]);
saveas(gcf,'Fig-cov-foursubfigs-first2last2-1960','epsc');
%savefig('Fig-cov-foursubfigs-first2last2-1960.fig')


