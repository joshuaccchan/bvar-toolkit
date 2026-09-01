clc; clear all;
addpath('./results_mat')
n = 20;
re_kappa = [];
load('fullsampleCS1.mat'); CS1 = Mean_Impact; re_kappa = [re_kappa, Mean_kappa];
load('fullsampleCS2.mat'); CS2 = Mean_Impact; re_kappa = [re_kappa, Mean_kappa];
load('fullsampleOI1.mat'); OI1 = Mean_Impact; re_kappa = [re_kappa, Mean_kappa];
load('fullsampleOI2.mat'); OI2 = Mean_Impact; re_kappa = [re_kappa, Mean_kappa];

mincolor = min([CS1(:);OI1(:)]); maxcolor = max([CS1(:);OI1(:)]); 

figure; subplot(1,2,1);
heatmap(CS1,[1:n],[1:n],[],'ColorMap', @cool,  ...
    'Colormap', 'money','MinColorValue', mincolor, 'MaxColorValue', maxcolor); title('CS');
subplot(1,2,2);
heatmap(OI1,[1:n],[1:n],[],'ColorMap', @cool, 'Colorbar', true, ...
    'Colormap', 'money','MinColorValue', mincolor, 'MaxColorValue', maxcolor); title('OI');
set(gcf,'position',[100,100,800,400]) 
%--- WARNING: need to hand-adjust the heatmap to make values of 0 pure white.
saveas(gcf,'hm-CS1OI1-start1960','epsc')


mincolor = min([CS2(:);OI2(:)]);
maxcolor = max([CS2(:);OI2(:)]);
figure; subplot(1,2,1);
heatmap(CS2,[1:n],[1:n],[],'ColorMap', @cool,  ...
    'Colormap', 'money','MinColorValue', mincolor, 'MaxColorValue', maxcolor); title('CS-2');
subplot(1,2,2);
heatmap(OI2,[1:n],[1:n],[],'ColorMap', @cool, 'Colorbar', true, ...
    'Colormap', 'money','MinColorValue', mincolor, 'MaxColorValue', maxcolor); title('OI-2');
set(gcf,'position',[100,100,800,400])
% saveas(gcf,'hm-CS2OI2-start1960','epsc')

re_kappa_Mean=re_kappa; 
csvwrite('Table2_kappa.csv',re_kappa_Mean);
