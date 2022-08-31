% 6/23/22 plot results of downsampling analysis on chr22

save_path = '/Volumes/broad_oconnor/trees/processed_results/main_results.txt';
T = readtable(save_path);
addpath('plotting_scripts')

% Filter by parameter settings
afr_rows = find(strcmp(T.pop,'AFR'));
[~,representatives] = unique(T.start(afr_rows));

rows = T.l1_pen == 0.1;
rows(afr_rows(representatives)) = true;

T = T(rows, :);

% Filter by blocksize
small_block = 59565357;
rows = T.start == small_block;
T = T(~rows, :);

% Median statistics for each population
[medians, ~, rows] = unique(T(:,6));
xdata = [];
for ii=1:size(medians,1)
    medians.mse(ii) = median(T.mse(rows == ii));
    medians.degree(ii) = median(T.degree(rows == ii));
    medians.initial_degree(ii) = median(T.initial_degree(rows == ii));
    medians.runtime(ii) = median(T.runtime(rows == ii));
    medians.maxruntime(ii) = max(T.runtime(rows == ii));
    medians.size(ii) = median(T.size(rows == ii));
    medians.count(ii) = sum(rows == ii);
    medians.meanruntime(ii) = mean(T.runtime(rows == ii));
    xdata = [xdata; ii* ones(sum(rows == ii),1)];
end
disp(medians);

colors = colororder();
color1 = colors(1,:);
outlier_color1 = 1 - (0.25 * (1 - color1));
color2 = colors(3,:);
outlier_color2 = 1 - (0.25 * (1 - color2));
figure
subplot(1,3,2)
boxchart(xdata,T.mse,'MarkerStyle', '.','boxfacecolor',color2,'markercolor',outlier_color2,'jitteroutliers',0);
% ymax = max(T.mse) * 1.01;
ymax = 0.004;
ylim([0 ymax])
ylabel('In-sample MSE')
set(gca,'XTick',1:max(xdata),'XTickLabel',...
    medians.pop,...
    'YTick',0:0.001:ymax,...
    'YTickLabel',arrayfun(@(x)sprintf('%.3f',x),0:0.001:ymax,'UniformOutput',0))
xticklabel_rotate([],45)


subplot(1,3,1)
boxchart(xdata,T.degree,'MarkerStyle', '.','boxfacecolor',color1,'markercolor',outlier_color1);
ymax = max(T.degree) * 1.01;
ylim([0 ymax])
ylabel('Edges per SNP')
set(gca,'XTick',1:max(xdata),'XTickLabel',...
    medians.pop,...
    'YTick',0:5:ymax)
xticklabel_rotate([],45)

save_path = '/Volumes/broad_oconnor/trees/processed_results/ukb_mse.txt';
T = readtable(save_path);

% Median statistics for each population
[medians, ~, rows] = unique(T(:,3));
xdata = [];
for ii=1:size(medians,1)
    medians.mse_PvR(ii) = median(T.mse_PvR(rows == ii));
    medians.mse_RvR(ii) = median(T.mse_RvR(rows == ii));
    xdata = [xdata; ii* ones(sum(rows == ii),1)];
end
disp(medians);

subplot(1,3,3)
boxchart(xdata,T.mse_PvR,'MarkerStyle', '.','boxfacecolor',color2,'markercolor',outlier_color2);
% ymax = max(T.mse_PvR) * 1.01;
ymax = 0.004;
ylim([0 ymax])
ylabel('pan-UKB MSE')
set(gca,'XTick',1:max(xdata),'XTickLabel',...
    medians.pop,...
    'YTick',0:0.001:ymax,...
    'YTickLabel',arrayfun(@(x)sprintf('%.3f',x),0:0.001:ymax,'UniformOutput',0))
xticklabel_rotate([],45)

