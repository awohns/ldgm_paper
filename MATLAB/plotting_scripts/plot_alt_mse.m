save_path = '/Volumes/broad_oconnor/trees/processed_results/altmse.txt';
T = readtable(save_path);

% Median statistics for each population
[medians, ~, rows] = unique(T(:,4));
xdata = [];
for ii=1:size(medians,1)
    medians.mse(ii) = median(T.mse(rows == ii));
    medians.altmse(ii) = median(T.mse(rows == ii));
    medians.size(ii) = median(T.size(rows == ii));
    medians.count(ii) = sum(rows == ii);
    xdata = [xdata; ii* ones(sum(rows == ii),1)];
end
disp(medians);

colors = colororder();
color1 = colors(1,:);
outlier_color1 = 1 - (0.25 * (1 - color1));
color2 = colors(3,:);
outlier_color2 = 1 - (0.25 * (1 - color2));

figure
boxchart(xdata,T.altmse,'MarkerStyle', '.','boxfacecolor',color2,'markercolor',outlier_color2,'jitteroutliers',0);
ymax = max(T.altmse)*1.01;
ylim([0 ymax])
ylabel('alternative MSE')
set(gca,'XTick',1:max(xdata),'XTickLabel',...
    medians.pop,...
    'YTick',0:1e-5:ymax)
