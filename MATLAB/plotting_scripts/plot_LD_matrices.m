clear;clc
addpath('plotting_functions')
filepath = '~/Dropbox/LD_graphical_models/data/';
filename = '1kg_ALL_chr22_22004675_23374984_MAF_0.01_RF_0.01_T_8';
X = load([filepath, filename, '.genos'])';
T = readtable([filepath, '1kg_nygc_trios_removed_All_pops_geno_ids_pops.csv']);

r2_threshold = 0.01;

superpops = table2cell(T(:,3));
[superpop_names, ~, which_pop] = unique(superpops);
% incl = which_pop == 4; %EUR

popns = [1 4 2 ];
included_SNPs = false(size(X,2), length(popns));

for ii=1:length(popns)
    % Allele frequency for each population
    incl = which_pop == popns(ii);
    af(:,ii) = mean(X(incl,:));
    
    % Precision matrix and its inverse for each population
    load([filepath, filename, '.path_distance=4.0.l1_pen=0.10.maf=0.01.',...
        superpop_names{popns(ii)},'.precisionMatrix.mat'], 'SNPs', ...
        'precisionEstimate', 'A_weighted');
    included_SNPs(SNPs,ii) = true;
    P{ii} = precisionEstimate;
    P_inv{ii} = inv(precisionEstimate);
end

% SNPs that are common + have precision matrix entry in every popn
present_all_popns = all(included_SNPs,2);

clear *degree
figure;
for ii=1:length(popns)
    incl = which_pop == popns(ii);
    subplot(2,5,ii)
    R =  corr(X(incl,present_all_popns));
    imagesc(triu(R));
    colormap(bluewhitered(256));
    caxis([-1 1])
    title(superpop_names(popns(ii)))
    set(gca,'XTick',[],'YTick',[])
    axis off
    box off
    
    R_degree(:,ii) = full(sum(R.^2 > r2_threshold));

    
    subplot(2,5,ii + 5)
    incl = present_all_popns(included_SNPs(:,ii));
    R = P_inv{ii}(incl, incl);
    imagesc(triu(R));
    colormap(bluewhitered(256));
    caxis([-1 1])
    %     title(superpop_names(popns(ii)))
    set(gca,'XTick',[],'YTick',[])
    axis off
    box off

    degree(:,ii) = full(sum(P{ii}(incl,incl) ~= 0));

end

nrows = size(degree,1);
xdata = ones(nrows,1) .* [1 3 5 2 4 6];
cdata = ones(nrows,1) .* [1 1 1 0 0 0];

data = [degree, R_degree];
figure;boxchart(xdata(:), data(:),'markerstyle','.' ,'GroupByColor', cdata(:))
%  ,'GroupByColor',xdata(:)
set(gca,'XTick',1.5: 2: 7.5, 'XTickLabel', superpop_names(popns))
legend('LD neighbors per SNP (r^2 > 0.01)', 'Precision matrix neighbors per SNP')
hold on
% line([0,9],[20,20],'color','black')
ylabel('Neighbors per SNP')
ylim([0 1000])
legend boxoff
% view([90 -90])
% figure;
% mymap =[1, 1, 1; 0.5, 0.5, 0.5; 1, 0, 0];
% A = A_weighted(incl, incl) > threshold;
% A = P{end}(incl,incl) - 0.5 * A;
% 
% subplot(2,1,1)
% imagesc(triu(A));
% colormap(mymap);  caxis([0 1])
% set(gca,'XTick',[],'YTick',[])
% 
% subplot(2,1,2)
% imagesc(triu(A(end-199:end,end-199:end)))
% colormap(mymap);  caxis([0 1])
% set(gca,'XTick',[],'YTick',[])


