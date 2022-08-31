save_path = '/Volumes/broad_oconnor/trees/processed_results/BLUP_UKBvs1kg.txt';
phenoTable_save_path = '/Volumes/broad_oconnor/trees/processed_results/BLUP_UKBvs1kg_pvals.txt';

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

figure
boxchart(xdata(:),ydata(:),'MarkerStyle', '.',...
    'GroupByColor',cdata(:));
set(gca,'XTick',1:4,'XTickLabel',pheno_names(phenoOrder))
ylabel('PGS r^2, 1kg LD vs UKB LD')
legend('BLUP-ldcov','BLUP-ldgm')
legend boxoff

hold on
for ii = 1:noPhenos
    text(ii,1,sprintf('p=%.0e',phenoT.pval_weights(phenoOrder(ii))))
end