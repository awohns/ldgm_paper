save_path = '/Volumes/broad_oconnor/luke/results/070122/file_size_table';

T = readtable(save_path);
xdata = T.blockSize;
ydata = [T.precisionSize T.correlationSize];
figure; hold on
for ii = 1:2
    scatter(xdata,ydata(:,ii))
end
xlabel('No. SNPs')
ylabel('File size')
legend('LDGM precision matrix','Correlation matrix')

[pops, ~, idx] = unique(T.pop);

for ii = 1:length(pops)
    precision_sum(ii) = sum(T.precisionSize(idx == ii))/1e9;    
    correlation_sum(ii) = sum(T.correlationSize(idx == ii))/1e9;
end

yy = [precision_sum' correlation_sum'];
figure
bar(yy)
set(gca,'XTickLabel',pops)

hold on

for ii=1:length(pops)
    str = sprintf('%.0fx',correlation_sum(ii)/precision_sum(ii));
    text(ii,correlation_sum(ii)+1,str);
end
box off
xdata = T.blockSize;
ydata = [T.precisionSize T.correlationSize];
figure; hold on
for ii = 1:2
    scatter(xdata,ydata(:,ii))
end
xlabel('No. SNPs')
ylabel('File size')
legend('LDGM precision matrix','Correlation matrix')

[pops, ~, idx] = unique(T.pop);

for ii = 1:length(pops)
    precision_sum(ii) = sum(T.precisionSize(idx == ii))/1e9;    
    correlation_sum(ii) = sum(T.correlationSize(idx == ii))/1e9;
end

yy = [precision_sum' correlation_sum'];
figure
bar(yy)
set(gca,'XTickLabel',pops)

hold on

for ii=1:length(pops)
    str = sprintf('%.0fx',correlation_sum(ii)/precision_sum(ii));
    text(ii,correlation_sum(ii)+1,str);
end
box off

%diffpcm = [population_names,(precision_sum)',(correlation_sum)'],

population_names = {'AFR';'AMR';'EAS';'EUR';'SAS'};
precision_size = (precision_sum)';
correlation_size = (correlation_sum)';


size_table = table(correlation_size,precision_size,...
    'RowNames',population_names)

fig = uifigure;
uit = uitable(fig,'Data',size_table);

