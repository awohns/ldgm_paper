% This script calculates BLUPx r2 for different sample sizes in different
% ancestry groups 
% 
% Use it to reproduce results from Figure 5d
% 
% Updated 8/14/22

clear;clc
jobnumber = str2num(getenv('SGE_TASK_ID'));
if isempty(jobnumber);jobnumber=1;end
rng(jobnumber);

addpath(genpath('/broad/oconnor/trees/ldgm/'))
ldgms_dir = '/broad/oconnor/trees/LDGMs_aug2022/';
populations = {'AFR','EUR','EAS','SAS'};
save_file = '/broad/oconnor/trees/results_081522/BLUPx_equalR2_081522.txt';

[P, snplists] = loadLDGMs(ldgms_dir,populations);
incl = ~any(cellfun(@isempty,P),2);
P = P(incl,:); snplists = snplists(incl);
[noBlocks, noPopns] = size(P);
disp(noBlocks)

AF = cell(noBlocks,noPopns);
for ii = 1:noBlocks
    [~,representatives] = unique(snplists{ii}.index,'stable');
    
    AF(ii,:) = mat2cell(table2array(snplists{ii}(representatives,populations)),...
        length(representatives),ones(1,noPopns));
end

crosspop_heritability = 0.1 * (eye(noPopns)*0.1 + ones(noPopns)*0.9);

sampleSize_array = [ones(1,noPopns)/4;
    0.7 0.1 0.1 0.1];
sampleSize_array = [1e4 * sampleSize_array;
    1e5 * sampleSize_array;
    1e6 * sampleSize_array];


r2_PGS = zeros(size(sampleSize_array,1), noPopns);
for sim = 1:size(sampleSize_array,1)
    [sumstats, whichIndices, prior_variance, beta_perallele] = ...
        simulateSumstats(sampleSize_array(sim,:), AF,...
        'precisionMatrices',P,...
        'heritability',crosspop_heritability);
    
    betaBLUPx = BLUPxldgm(P, whichIndices, sumstats, prior_variance, sampleSize_array(sim,:));
    
    r2_PGS(sim,:) = r2PGS(beta_perallele,betaBLUPx,P,whichIndices,AF);
    
end

fs = [repmat('%f ',1,noPopns),'\n'];
file = fopen(save_file,'a');
fprintf(file,fs,r2_PGS);
fclose(file);









