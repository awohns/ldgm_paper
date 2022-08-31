clear;clc
save_path = '/Volumes/broad_oconnor/trees/processed_results/matrix_operations_runtime.txt';
T = readtable(save_path);


labels = T.description;

time = [T.time_ldcov T.time_ldgm]/60^2;

incl = [8 7 5 4 3 6 2 1];
figure;
barh(time(incl,:))
set(gca,'YTickLabel',labels(incl))
xmax = max(time(incl,1))*1.05;
xlim([0 xmax])
set(gca,'XTick',0:3:xmax)
xlabel('Runtime (hours)')
plots = get(gca,'Children');
legend(plots,'LDGM precision matrix','Correlation matrix')
legend boxoff

for ii=1:length(incl)
    str = sprintf('%.0fx',(time(incl(ii),1)/time(incl(ii),2)));
    text(max(time(incl(ii),:))+0.1,ii,str);
end
box off

% 
% figure
% subplot(1,2,1)
% incl = [6 2 1];
% 
% barh(time(incl,:))
% set(gca,'YTickLabel',labels(incl))
% xlim([0 max(time(incl,2))*1.05])
% xlabel('Runtime (minutes)')
% 
% for ii=1:length(incl)
%     str = sprintf('%.0fx',(time(incl(ii),2)/time(incl(ii),1)));
%     text(max(time(incl(ii),:))+1,ii,str);
% end
% box off
% 
% subplot(1,2,2)
% incl = [10 5 4 3];
% 
% barh(time(incl,:))
% set(gca,'YTickLabel',labels(incl))
% xlim([0 max(time(incl,2))*1.05])
% xlabel('Runtime (minutes)')
% plots = get(gca,'Children');
% legend(plots,'LDGM precision matrix','Correlation matrix')
% 
% for ii=1:length(incl)
%     str = sprintf('%.0fx',(time(incl(ii),2)/time(incl(ii),1)));
%     text(max(time(incl(ii),:))+1,ii,str);
% end
% legend boxoff
% box off





