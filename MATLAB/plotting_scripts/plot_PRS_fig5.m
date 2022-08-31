clear;clc
h2 = 0.1;
runtime_path = '/Volumes/broad_oconnor/trees/processed_results/BLUPx_runtime.txt';
equal_r2_path = '/Volumes/broad_oconnor/trees/processed_results/BLUPx_equalR2_081522.txt';
asym_r2_path = '/Volumes/broad_oconnor/trees/processed_results/asym_nonport_081622.txt';
save_path = '/Volumes/broad_oconnor/trees/processed_results/BLUP_UKBvs1kg.txt';
phenoTable_save_path = '/Volumes/broad_oconnor/trees/processed_results/BLUP_UKBvs1kg_pvals.txt';
titles = {'a In-sample vs reference LD', 'b BLUPx Runtime', 'c Train in one population', 'd Train in four populations'};

figure
%% Fig 5a
subplot(2,2,1)

T = readtable(save_path);
phenoT = readtable(phenoTable_save_path);
[phenos, ~, idx] = unique(T.phenotype);
pheno_names = {'BMI','CVD','Height','T2D'};
noPhenos = length(phenos);
noBlocks = height(T) / noPhenos;
phenoOrder = [3 1 2 4];
xdata = zeros(2 * noBlocks, noPhenos); ydata = xdata; 
cdata = [zeros(noBlocks, noPhenos); ones(noBlocks, noPhenos);];
for ii = 1:noPhenos
    xdata(:,ii) = ii;
    ydata(:,ii) = [T.ldcov_weights_r2(idx == phenoOrder(ii)); T.ldgm_weights_r2(idx == phenoOrder(ii))];
end

boxchart(xdata(:),ydata(:),'MarkerStyle', '.',...
    'GroupByColor',cdata(:));
set(gca,'XTick',1:4,'XTickLabel',pheno_names(phenoOrder))
ylabel('r^2 (\beta_{UKB} with \beta_{1kg})')
legend('BLUP-ldcov','BLUP-ldgm')
legend boxoff

hold on
for ii = 1:noPhenos
    text(ii,1,sprintf('p=%.0e',phenoT.pval_weights(phenoOrder(ii))))
end

title(titles{1})

%% Fig 5b
Time = readtable(runtime_path);
disp(Time)
subplot(2,2,2)
plot2 = barh(1:5, [Time.ldcov_mean Time.ldgm_mean]/3600);
set(gca,'XTick',0:6:max(Time.ldcov_mean)/3600)
xlabel('Runtime (hours)')
legend(plot2(1:2),{  'BLUPx-ldcov', 'BLUPx-ldgm',})
legend boxoff
set(gca,'YTIck',1:5,'YTickLabel',5:-1:1)
ylabel('No. populations')
xlim([0 1.01*max(Time.ldcov_mean)/3600])
for ii = 1:5
    text(Time.ldcov_mean(ii)/3600,ii,sprintf('%dx',round(Time.ldcov_mean(ii)/Time.ldgm_mean(ii))));
end
title(titles{2})

%% Fig 5c
canonical_cmap = [253, 232, 11; 107, 174, 213; 57, 161, 83; 123, 110, 176]/256;

pop1 = 'AFR';
pop2 = 'EUR';
popns = {pop1, pop2};

data = readtable(asym_r2_path);
nn = 1e5;
h2 = 0.1;

for p1 = 1:2
    for p2 = 1:2
        plotvals(p1,p2) = h2 * data.r2(strcmp(data.training,popns{p1}) & ...
            strcmp(data.testing,popns{p2}) & data.sample_size == nn);
    end
end

subplot(2,3,4);
% set(gcf,'Position',[418 321 378 420]);

h = bar(plotvals,'grouped');

hold on;
ngroups = size(plotvals, 1);
nbars = size(plotvals, 2);

cat_label = categorical({'Train in AFR', 'Train in EUR'}); %The Group Label
set(gca,'xticklabel',cat_label,'XTickLabelRotation',0)
ylabel('r^2_{PGS}');

ytickmat = 0:0.02:0.09;
set(gca,'YTick',0:0.02:0.1,'YTickLabel',[mat2cell(ytickmat,1,length(ytickmat)), {'h^2'}])
hold on
line([0.5 2.5],[h2 h2],'linestyle','--','color','black')
ylim([0 1.1*h2])
xlim([0.5 2.5])
for popn = 1:2
    h(popn).FaceColor = canonical_cmap(popn,:);
end
% legend(h,{['Target=',pop1],['Target=',pop2]},'Location','southoutside')
% legend boxoff
title(titles{3})

%% Fig 5d

popns = {'AFR','EUR','EAS','SAS'};

nn_label = {'25% each',...
    '70% AFR, 10% others',...
    };
nn_array_order = [5:6]; % Make sure this corresponds with the order
                      % in the simulation script
h2 = 0.1;
noPops = length(popns);
noNs = length(nn_array_order);

data = readtable(equal_r2_path);

for i1 = 1:noNs
    for i2 = 1:noPops
        plotvals(i1,i2) =  h2 * data.r2mean(data.sim == nn_array_order(i1) & strcmp(data.popn,popns(i2)));
        err(i1,i2) =  h2 * data.r2sd(data.sim == nn_array_order(i1) & strcmp(data.popn,popns(i2)));
    end
end

subplot(2,2,4);

b = bar(plotvals);
for popn = 1:length(popns)
    b(popn).FaceColor = canonical_cmap(popn,:);
end
legend([b],cellfun(@(p){['Predict in ',p]},popns),'Orientation','horizontal','Location','southoutside','Autoupdate','off');
legend boxoff

% Adding error bars
hold on;
ngroups = size(plotvals, 1);
nbars = size(plotvals, 2);
% Calculating the width for each bar group
% groupwidth = min(0.8, nbars/(nbars + 1.5));
% for i = 1:nbars
%     x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
%     errorbar(x, plotvals(:,i), err(:,i), 'k.');
% end
hold off

ytickmat = 0:0.02:0.09;
set(gca,'YTick',0:0.02:0.1,'YTickLabel',[mat2cell(ytickmat,1,length(ytickmat)), {'h^2'}])
line([0.5 2.5],[h2,h2],'linestyle','--','color','black')
ylim([0 1.1*h2])

set(gca,'xticklabel',nn_label,'XTickLabelRotation',0)
ylabel('r^2_{PGS}')
% xlabel('Training Population Composition');
title(titles{4})