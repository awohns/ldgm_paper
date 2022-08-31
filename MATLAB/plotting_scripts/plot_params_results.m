% 6/23/22 plot results of analysis varying precision inference parameters on chr22

save_path = '/Volumes/broad_oconnor/trees/processed_results/params_chr22.txt';
T = readtable(save_path);

[medians, ~, rows] = unique(T(:,4:5));
for ii=1:size(medians,1)
    medians.mse(ii) = median(T.mse(rows == ii));
    medians.degree(ii) = median(T.degree(rows == ii));
    medians.initial_degree(ii) = median(T.initial_degree(rows == ii));
    medians.runtime(ii) = median(T.runtime(rows == ii));
    medians.maxruntime(ii) = max(T.runtime(rows == ii));
    medians.count(ii) = sum(rows == ii);
end
disp(medians);

mse_data = [T(T.l1_pen == 0.1 & T.distance_threshold==3,:).mse;...
    T(T.l1_pen == 0.1 & T.distance_threshold==4,:).mse;...
    T(T.l1_pen == 0.1 & T.distance_threshold==5,:).mse;...
    T(T.l1_pen == 0.1 & T.distance_threshold==6,:).mse];

degree_data = [T(T.l1_pen == 0.1 & T.distance_threshold==3,:).degree;...
    T(T.l1_pen == 0.1 & T.distance_threshold==4,:).degree;...
    T(T.l1_pen == 0.1 & T.distance_threshold==5,:).degree;...
    T(T.l1_pen == 0.1 & T.distance_threshold==6,:).degree];

xdata = [];
for ii=3:6; 
    njobs(ii) = sum(T.l1_pen == 0.1 & T.distance_threshold==ii); 
    xdata = [xdata; (ii-2) * ones(njobs(ii),1)];
end


figure; 
subplot(2,2,1)
boxchart(xdata(:),mse_data(:),'MarkerStyle', '.');
ymax = max(mse_data(:)) * 1.01;
ylim([0 ymax])
ylabel('Mean sq error')
set(gca,'XTick',1:4,'XTickLabel',...
    3:6,...
    'YTick',0:0.001:ymax,'YTickLabel',arrayfun(@(x)sprintf('%.3f',x),0:0.001:ymax,'UniformOutput',0))
xlabel('Distance threshold')
title('a MSE by path distance threshold')

subplot(2,2,2)
boxchart(xdata(:),degree_data(:),'MarkerStyle', '.');
ymax = max(degree_data(:)) * 1.01;
ylim([0 ymax])
ylabel('Edges per SNP')
set(gca,'XTick',1:4,'XTickLabel',...
    3:6,...
    'YTick',0:10:ymax,'YTickLabel',0:10:ymax)
xlabel('Distance threshold')
title('b Density by path distance threshold')

mse_data = [T(T.l1_pen == 0 & T.distance_threshold==4,:).mse;...
    T(T.l1_pen == 0.05 & T.distance_threshold==4,:).mse;...
    T(T.l1_pen == 0.1 & T.distance_threshold==4,:).mse;...
    T(T.l1_pen == 0.2 & T.distance_threshold==4,:).mse];

degree_data = [T(T.l1_pen == 0 & T.distance_threshold==4,:).degree;...
    T(T.l1_pen == 0.05 & T.distance_threshold==4,:).degree;...
    T(T.l1_pen == 0.1 & T.distance_threshold==4,:).degree;...
    T(T.l1_pen == 0.2 & T.distance_threshold==4,:).degree];

xdata = [ones(sum(T.l1_pen == 0 & T.distance_threshold==4),1);...
    2 * ones(sum(T.l1_pen == 0.05 & T.distance_threshold==4),1);...
    3 * ones(sum(T.l1_pen == 0.1 & T.distance_threshold==4),1);...
    4 * ones(sum(T.l1_pen == 0.2 & T.distance_threshold==4),1)];

subplot(2,2,3)
boxchart(xdata(:),mse_data(:),'MarkerStyle', '.');
ymax = max(mse_data(:)) * 1.01;
ylim([0 ymax])
ylabel('Mean sq error')
set(gca,'XTick',1:4,'XTickLabel',...
    [0 0.05 0.1 0.2],...
    'YTick',0:0.001:ymax,'YTickLabel',arrayfun(@(x)sprintf('%.3f',x),0:0.001:ymax,'UniformOutput',0))
xlabel('l1 penalty')
title('c MSE by l1 penalty')

subplot(2,2,4)
boxchart(xdata(:),degree_data(:),'MarkerStyle', '.');
ymax = max(degree_data(:)) * 1.01;
ylim([0 ymax])
ylabel('Edges per SNP')
set(gca,'XTick',1:4,'XTickLabel',...
    [0 0.05 0.1 0.2],...
    'YTick',0:10:ymax,'YTickLabel',0:10:ymax)
xlabel('l1 penalty')
title('d Density by l1 penalty')

























