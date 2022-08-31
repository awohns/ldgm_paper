% This script runs and times BLUPx-ldgm and BLUPx-ldcov for different 
% numbers of ancestry groups 
% 
% Use it to reproduce results from Figure 5b
% 
% Updated 8/24/22

clear;clc
jobnumber = str2num(getenv('SGE_TASK_ID'));
if isempty(jobnumber);jobnumber=1;end
rng(jobnumber);

addpath(genpath('/broad/oconnor/trees/ldgm/'))
ldgms_dir = '/broad/oconnor/trees/LDGMs_aug2022/';
populations = {'AFR','EUR','EAS','SAS'};
save_file = '/broad/oconnor/trees/results_081522/BLUPx_runtime.txt';

[P, snplists] = loadLDGMs(ldgms_dir,populations,jobnumber);
incl = ~any(cellfun(@isempty,P),2);
P = P(incl,:); snplists = snplists(incl);
[noBlocks, noPopns] = size(P);

AF = cell(noBlocks,noPopns);
for ii = 1:noBlocks
    [~,representatives] = unique(snplists{ii}.index,'stable');
    
    AF(ii,:) = mat2cell(table2array(snplists{ii}(representatives,populations)),...
        length(representatives),ones(1,noPopns));
end

crosspop_heritability = 0.1 * (eye(noPopns)*0.1 + ones(noPopns)*0.9);

sampleSize_array = 1e4 * triu(ones(noPopns))';

for ii = 1:noBlocks
    for jj = 1:noPopns
        whichIndices{ii,jj} = diag(P{ii,jj} ~= 0);
        allTrue{ii,jj} = true(sum(whichIndices{ii,jj}),1);
        R{ii,jj} = inv(P{ii,jj}(whichIndices{ii,jj},whichIndices{ii,jj}));
    end
end

runtime = zeros(size(sampleSize_array,1), 2);
for sim = 1:size(sampleSize_array,1)
    include_popn = sampleSize_array(sim,:) ~=0;
    [sumstats, ~, prior_variance, beta_perallele] = ...
        simulateSumstats(sampleSize_array(sim,include_popn), AF(:,include_popn),...
        'precisionMatrices',P(:,include_popn),...
        'heritability',crosspop_heritability(include_popn,include_popn));
    
    tic;
    BLUPxldgm(P(:,include_popn), whichIndices(:,include_popn), sumstats, prior_variance, sampleSize_array(sim,include_popn));
    runtime(sim,1) = toc;

    tic;
    BLUPxldcov(R(:,include_popn), allTrue(:,include_popn), sumstats, prior_variance, sampleSize_array(sim,include_popn));
    runtime(sim,2) = toc
%     r2_PGS(sim,:) = r2PGS(beta_perallele,betaBLUPx,P,whichIndices,AF);
    
end

fs = [repmat('%f ',1,numel(runtime)),'\n'];
file = fopen(save_file,'a');
fprintf(file,fs,runtime);
fclose(file);









