% 6/23/22 plot results of downsampling analysis on chr22

save_path = '/Volumes/broad_oconnor/trees/processed_results/downsample_chr22.txt';
T = readtable(save_path);

data = [T.mse T.mse_crossval T.mse_RvsR];
xdata = repmat(1:3, size(T,1), 1);

figure
boxchart(xdata(:),data(:),'MarkerStyle', '.');
ymax = max(data(:)) * 1.01;
ylim([0 ymax])
ylabel('Mean sq error')
set(gca,'XTick',1:3,'XTickLabel',...
    {'P vs R in-sample','P vs R cross-val','R vs R cross-val'},...
    'YTick',0:0.001:ymax,'YTickLabel',arrayfun(@(x)sprintf('%.3f',x),0:0.001:ymax,'UniformOutput',0))
xticklabel_rotate([],45)
title('MSE with 50% samples held out')
