clear;clc

save_path = '/Volumes/broad_oconnor/trees/processed_results/MAFbin_MSE.txt';
T = readtable(save_path);

figure
subplot(1,2,1)
boxchart([T.mse_1 T.mse_2 T.mse_3],'markerstyle','.')
ymax = 0.002;
ylim([0 ymax])
set(gca,'XTickLabel',...
    {'0.01-0.05','0.05-0.2','0.2-0.5'},...
    'YTick',0:0.001:ymax,...
    'YTickLabel',arrayfun(@(x)sprintf('%.3f',x),0:0.001:ymax,'UniformOutput',0))
xlabel('MAF range')
ylabel('MSE')
title('a. Accuracy by MAF bin')
% xticklabel_rotate([],45)

subplot(1,2,2)
boxchart([T.degree_1 T.degree_2 T.degree_3],'markerstyle','.')
ymax = max(T.degree_3)*1.01;
ylim([0 ymax])
set(gca,'XTickLabel',...
    {'0.01-0.05','0.05-0.2','0.2-0.5'},...
    'YTick',0:5:ymax)
xlabel('MAF range')
ylabel('Edges per SNP')
title('b. Density by MAF bin')
