
save_path = '/Volumes/broad_oconnor/trees/processed_results/';
save_path_ldgm = [save_path,'existing_approaches_ldgm.txt'];
save_path_lowrank = [save_path,'existing_approaches_lowrank.txt'];

Tldgm = readtable(save_path_ldgm);
Tlowrank = readtable(save_path_lowrank);
noFiles = height(Tldgm);
kk = unique(Tlowrank.rank,'stable');
noKValues = numel(kk);
assert(height(Tlowrank) == noFiles * noKValues);

figure;subplot(3,1,1)
data = [reshape(Tlowrank.mse,noFiles,noKValues) Tldgm.mse];
xgroupdata = repmat(1:size(data,2),size(data,1),1);
cgroupdata = [ones(noFiles,noKValues) 2*ones(noFiles,1)];
b = boxchart(xgroupdata(:),data(:),'GroupByColor',cgroupdata(:),'MarkerStyle','.');
ylabel('MSE')
set(gca,'XTick',[(1:size(data,2)-1)-0.25, size(data,2)+0.25],'XTickLabel',[{'k=5'},num2cell(kk(2:end)'),'LDGM'])
set(gca,'YTick',10.^-(8:-2:2))
set(gca,'YScale','log')
hold on;title('a Mean-squared error of R_k vs R_{LDGM}')

subplot(3,1,2)
data = [reshape(Tlowrank.altmse,noFiles,noKValues) Tldgm.altmse];
xgroupdata = repmat(1:size(data,2),size(data,1),1);
cgroupdata = [ones(noFiles,noKValues) 2*ones(noFiles,1)];
b = boxchart(xgroupdata(:),data(:),'GroupByColor',cgroupdata(:),'MarkerStyle','.');
ylabel('alt MSE ratio')
set(gca,'XTick',[(1:size(data,2)-1)-0.25, size(data,2)+0.25],'XTickLabel',[{'k=5'},num2cell(kk(2:end)'),'LDGM'])
set(gca,'YTick',10.^-(8:-2:-2))
set(gca,'YScale','log')
hold on;title('b altMSE ratio of P_k vs P_{LDGM}')

subplot(3,1,3)
data = [reshape(Tlowrank.mse_with_ldgm,noFiles,noKValues)];
xgroupdata = repmat(1:size(data,2),size(data,1),1);
b = boxchart(xgroupdata(:),data(:),'MarkerStyle','.');
ylabel('MSE')
set(gca,'XTick',[(1:size(data,2))],'XTickLabel',[{'k=5'},num2cell(kk(2:end)')])
set(gca,'YScale','log')
% xticklabel_rotate([],45);
hold on;title('c MSE between R_k and R_{ldgm}')


