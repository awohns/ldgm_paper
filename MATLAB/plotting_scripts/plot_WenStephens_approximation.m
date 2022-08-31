
save_path = '/Volumes/broad_oconnor/trees/processed_results/';
save_path_ldgm = [save_path,'existing_approaches_ldgm.txt'];
save_path_WS = [save_path,'existing_approaches_WenStephens.txt'];

Tldgm = readtable(save_path_ldgm);
TWS = readtable(save_path_WS);
noFiles = height(Tldgm);
ldgm_median_mse = median(Tldgm.mse);
kk = unique(TWS.samplesize_param,'stable');
noKValues = numel(kk);
assert(height(TWS) == noFiles * noKValues);

figure;subplot(3,1,1)
data = [reshape(TWS.mse_with_ldgm,noFiles,noKValues)];
xgroupdata = repmat(1:size(data,2),size(data,1),1);

hold on
b = boxchart(xgroupdata(:),data(:),'MarkerStyle','.');
ylabel('MSE')
set(gca,'XTick',[(1:size(data,2))],'XTickLabel',round(kk))
set(gca,'YTick',10.^-(6:-1:2))
xlabel('Sample size parameter')
line([0 7],[1 1]*ldgm_median_mse,'linestyle','--','color','black')
% xticklabel_rotate([],45);
set(gca,'YScale','log')
hold on;title('a Mean-squared error between R_{WS} and R_{LDGM}')

subplot(3,1,2)
data = [reshape(TWS.mse,noFiles,noKValues)];
xgroupdata = repmat(1:size(data,2),size(data,1),1);
b = boxchart(xgroupdata(:),data(:),'MarkerStyle','.');
ylabel('MSE')
set(gca,'XTick',[(1:size(data,2))],'XTickLabel',round(kk))
set(gca,'YTick',10.^-(6:-2:2))
xlabel('Sample size parameter')
set(gca,'YScale','log')
hold on;title('b Mean-squared error between R_{WS} and R')

subplot(3,1,3)
data = [reshape(TWS.avgDegree,noFiles,noKValues)];
xgroupdata = repmat(1:size(data,2),size(data,1),1);
b = boxchart(xgroupdata(:),data(:),'MarkerStyle','.');
ylabel('Average degree')
set(gca,'XTick',[(1:size(data,2))],'XTickLabel',round(kk))
xlabel('Sample size parameter')
set(gca,'YScale','log')
hold on;title('c R_{WS} nonzero entries per SNP')


