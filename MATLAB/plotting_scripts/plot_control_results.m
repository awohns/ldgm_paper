% 6/23/22 plot results of downsampling analysis on chr22

save_path = '/Volumes/broad_oconnor/trees/processed_results/control_chr22.txt';
T = readtable(save_path);

[medians, ~, rows] = unique(T(:,[6 4]));
for ii=1:size(medians,1)
    medians.mse(ii) = median(T.mse(rows == ii));
    medians.degree(ii) = median(T.degree(rows == ii));
    medians.initial_degree(ii) = median(T.initial_degree(rows == ii));
    medians.runtime(ii) = median(T.runtime(rows == ii));
    medians.maxruntime(ii) = max(T.runtime(rows == ii));
    medians.count(ii) = sum(rows == ii);
end
disp(medians);


mse_data = [T(strcmp(T.control,'') & T.distance_threshold==4,:).mse;...
    T(strcmp(T.control,'banded_control') & T.distance_threshold==4,:).mse;...
    T(strcmp(T.control,'r2_control') & T.distance_threshold==4,:).mse;...
    T(strcmp(T.control,'banded_control') & T.distance_threshold==8,:).mse];

degree_data = [T(strcmp(T.control,'') & T.distance_threshold==4,:).degree;...
    T(strcmp(T.control,'banded_control') & T.distance_threshold==4,:).degree;...
    T(strcmp(T.control,'r2_control') & T.distance_threshold==4,:).degree;...
    T(strcmp(T.control,'banded_control') & T.distance_threshold==8,:).degree];

xdata = [ones(size(T(strcmp(T.control,'') & T.distance_threshold==4,:),1),1);...
    2 * ones(size(T(strcmp(T.control,'banded_control') & T.distance_threshold==4,:),1),1);...
    3 * ones(size(T(strcmp(T.control,'r2_control') & T.distance_threshold==4,:),1),1);...
    4 * ones(size(T(strcmp(T.control,'banded_control') & T.distance_threshold==8,:),1),1)];

figure; 
subplot(1,2,1)
boxchart(xdata(:),mse_data(:),'MarkerStyle', '.');
ymax = max(mse_data(:)) * 1.01;
ylim([0 ymax])
ylabel('Mean sq error')
set(gca,'XTick',1:4,'XTickLabel',...
    {'Inferred LDGM','Banded LDGM','r^2 threshold LDGM','Wide banded LDGM'},...
    'YTick',0:0.005:ymax,'YTickLabel',arrayfun(@(x)sprintf('%.3f',x),0:0.005:ymax,'UniformOutput',0))
xticklabel_rotate([],45)
title('a. MSE vs control LDGMs')

subplot(1,2,2)
boxchart(xdata(:),degree_data(:),'MarkerStyle', '.');
ymax = max(degree_data(:)) * 1.01;
ylim([0 ymax])
ylabel('Edges per SNP')
set(gca,'XTick',1:4,'XTickLabel',...
    {'Inferred LDGM','Banded LDGM','r^2 threshold LDGM','Wide banded LDGM'},...
    'YTick',0:10:ymax,'YTickLabel',0:10:ymax)
xticklabel_rotate([],45)
title('b. Density vs control LDGMs')








