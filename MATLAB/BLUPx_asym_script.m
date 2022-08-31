% This script calculates BLUP r2 training in one population and testing in
% another, with different sample sizes 
% 
% Use it to reproduce results of figure 5c
% 
% Updated 8/14/22

clear;clc
jobnumber = str2num(getenv('SGE_TASK_ID'));
if isempty(jobnumber);jobnumber=1;end
rng(jobnumber);

addpath(genpath('/broad/oconnor/trees/ldgm/'))
ldgms_dir = '/broad/oconnor/trees/LDGMs_aug2022/';
populations = {'AFR','EUR','EAS','SAS'};
save_file = '/broad/oconnor/trees/results_081522/BLUP_asym_081522.txt';

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

sampleSize_array = [1e4 1e5 1e6]' .* ones(1,noPopns);

r2_PGS = zeros(noPopns,noPopns,size(sampleSize_array,1));
for sim = 1:size(sampleSize_array,1)
    [sumstats, whichIndices, prior_variance, beta_perallele] = ...
        simulateSumstats(sampleSize_array(sim,:), AF,...
        'precisionMatrices',P,...
        'heritability',crosspop_heritability);
    
    % Train in each population individually
    for popn = 1:noPopns
        betaBLUPx = BLUPxldgm(P(:,popn), whichIndices(:,popn),...
            sumstats(:,popn), prior_variance(popn,popn), sampleSize_array(sim,popn));
        
        r2_PGS(popn,:,sim) = r2PGS(beta_perallele,betaBLUPx,P,whichIndices,AF);
    end
end

r2_PGS = reshape(r2_PGS,noPopns^2,size(sampleSize_array,1)); 

fs = [repmat('%f ',1,size(sampleSize_array,1)),'\n'];
file = fopen(save_file,'a');
fprintf(file,fs,r2_PGS);
fclose(file);











