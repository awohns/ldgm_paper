load('../../data/inf_blockError.mat')
inferred = blockError;
load('../../data/sim_blockError.mat')
simulated = blockError;

figure
colors = colororder();
color1 = colors(1,:);
outlier_color1 = 1 - (0.25 * (1 - color1));
color2 = colors(3,:);
outlier_color2 = 1 - (0.25 * (1 - color2));

boxchart([ones(length(simulated),1); 2*ones(length(simulated),1)],[simulated'; inferred'],'MarkerStyle', '.','boxfacecolor',color2,'markercolor',outlier_color2)
set(gca,'XTick',1:2,'XTickLabel',{'True TS','Inferred TS'})
ymax = 0.004;
ylim([0 ymax])
set(gca,'YTick',0:.001:ymax,...
    'YTickLabel',arrayfun(@(x)sprintf('%.3f',x),0:0.001:ymax,'UniformOutput',0))
xticklabel_rotate([],45)
ylabel('In-sample MSE')
