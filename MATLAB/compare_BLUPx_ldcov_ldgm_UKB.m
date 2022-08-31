% Run BLUP-ldgm and BLUP-ldcov on UKB sumstats with UKB or 1kg LD 
% 
% Use it to reproduce results of Figure 5a
% 
% Updated 8/23/22

clear;clc
ii = 1;

addpath(genpath('/broad/oconnor/trees/ldgm/MATLAB'))
sumstats_fp = '/broad/oconnor/sumstats/imputed/';
precision_fp ='/broad/oconnor/trees/LDGMs_aug2022/';
genos_1kg_fp = '/broad/oconnor/trees/genos/';
sample_ids_fp = '/broad/oconnor/trees/LDGMs_nygc/1kg_nygc_trios_removed_All_pops_geno_ids_pops.csv';
ldgm_suffix = '.EUR';
UKB_cov_fp = '/broad/oconnor/trees/UKB/cov/';
% sumstats_filename = 'cov_EDU_YEARS.sumstats';
sumstats_filename_array = {'disease_T2D.sumstats',...
    'disease_CARDIOVASCULAR.sumstats',...
    'body_HEIGHTz.sumstats',...
    'body_BMIz.sumstats'};
pheno_name_array = {'Type II diabetes','Cardiovascular disease','Height','BMI'};
h2_array = [.073, .155, .570, .303]; % Loh et al 2018 Supplementary Table 3
N_array = [440 450 650 500] * 1e3;
save_path = '/broad/oconnor/trees/processed_results/BLUP_UKBvs1kg.txt';
phenoTable_save_path = '/broad/oconnor/trees/processed_results/BLUP_UKBvs1kg_pvals.txt';
chr = 21:22;

file_name_parser = @(s)textscan(s,...
    'ukb_chr%d_%d_%d.%*s');

ids_table = readtable(sample_ids_fp);
EUR_samples = strcmp(ids_table.id_superpops,'EUR');

files = [];
for cc = chr
    files = [files;dir([UKB_cov_fp, 'ukb_chr',num2str(cc), '_*.snplist'])];
end
noBlocks = length(files);

T = table('size',[noBlocks * length(N_array),0]);
phenoTable = table('size',[length(N_array),0]);

for pheno = 1:length(N_array)
    nn = N_array(pheno);
    h2 = h2_array(pheno);
    sumstats_filename = sumstats_filename_array{pheno};
    
    % Load sumstats
    sumstats = readtable([sumstats_fp,sumstats_filename],'filetype','text');
    total_heterozygosity = sum(2 * sumstats.EAF .* (1 - sumstats.EAF));
    betaCov = h2 ./ total_heterozygosity;
    
    sumstats = sumstats(ismember(sumstats.CHR, chr), :);
    
    % column of sumstats named 'EAF' needs to be renamed 'AF'
    AF_col = strcmp(sumstats.Properties.VariableNames,'EAF');
    sumstats.Properties.VariableNames{AF_col} = 'AF';

    for blocknumber = 1:noBlocks
        fields = file_name_parser(files(blocknumber).name);
        T.chr(blocknumber + noBlocks * (pheno - 1)) = fields{1};
        T.start(blocknumber + noBlocks * (pheno - 1)) = fields{2};
        T.stop(blocknumber + noBlocks * (pheno - 1)) = fields{3};
        T.phenotype(blocknumber + noBlocks * (pheno - 1)) = pheno_name_array(pheno);
        
        % Ridge parameter for LD correlation matrix
        lambda = 0;

        % UKB LD covariance matrix files
        filename = files(blocknumber).name(1:end-length('.snplist'));
        
        % UKB LD covariance data
        UKB_cov_snplist = readtable([UKB_cov_fp, filename, '.snplist'],'FileType','text');
        R_UKB = readmatrix([UKB_cov_fp, filename, '.correlation'],'FileType','text');
        
        % Load precision matrix that matches the file name
        [P,SNPlist] = loadLDGMs([precision_fp,'*',filename(5:end)],'EUR');
        assert(numel(P)==1)
        assert(all(size(P{1})==numel(unique(SNPlist{1}.index))))
        [~,representatives] = unique(SNPlist{1}.index);

        % calculate 1kg correlation matrix
        genos_1kg = readmatrix([genos_1kg_fp, '1kg_', filename(5:end), '.genos'],'FileType','text');
        assert(size(genos_1kg,1) == height(SNPlist{1}))
        R_1kg = corr(genos_1kg(representatives,EUR_samples)');
        
        % Merge sumstats with LD matrix SNPs
        [~, isumstats, iUKBSNPs ] = ...
            intersect(sumstats.SNP, UKB_cov_snplist.rsid, 'stable');
        mergedSumstats = sumstats(isumstats, :);
        R_UKB = R_UKB(iUKBSNPs, iUKBSNPs);

        % Merge SNPs between LDGM + sumstats
        mergedSumstats.row = (1:height(mergedSumstats))';
        [whichIndices,mergedSumstats,idx] = ...
            mergesnplists(SNPlist,mergedSumstats,P);
        whichIndices{1} = unfind(whichIndices{1},length(representatives));
        R_UKB = R_UKB(mergedSumstats{1}.row, mergedSumstats{1}.row);

        % Calculate MSE between each 1kg matrix and the UKB matrix
        nzd = diag(P{1})~=0;
        P_inv = inv(P{1}(nzd,nzd));
        P_inv = P_inv(whichIndices{1}(nzd),whichIndices{1}(nzd));
        T.ldcov_mse(blocknumber + noBlocks * (pheno - 1)) = mean((R_1kg(whichIndices{1},whichIndices{1}) - R_UKB).^2,'all');
        T.ldgm_mse(blocknumber + noBlocks * (pheno - 1)) = mean((P_inv - R_UKB).^2,'all');

        % Run BLUPx-ldgm
        tic;
        [~, betaBLUPxLDGM] =...
            BLUPxldgm(P, whichIndices, mergedSumstats, betaCov);
        timeLDGM = toc;
        
        % Run BLUPx-ldcov with 1kg LD data
        tic;
        [~, betaBLUPxLDcov] =...
            BLUPxldcov({R_1kg}, whichIndices, mergedSumstats, betaCov);
        timeLDcov = toc;
        
        % Run BLUPx-ldcov with UKB LD data
        tic;
        [~, betaBLUPxUKBLD] =...
            BLUPxldcov({R_UKB}, {(1:length(R_UKB))'}, mergedSumstats, betaCov);
        betaBLUPxUKBLD = cellfun(@assignto,betaBLUPxUKBLD,whichIndices,...
            {length(betaBLUPxLDcov{1})},'UniformOutput',false);
        timeLDcov = toc;
        
        % correlation between PGS weights for UKB LD covariance vs 1kg LD
        % covariance matrices
        T.ldcov_r2(blocknumber + noBlocks * (pheno - 1)) = corr(betaBLUPxLDcov{1},betaBLUPxUKBLD{1})^2;
        
        % correlation between PGS weights for UKB LD covariance vs 1kg LDGM
        % precision matrices
        T.ldgm_r2(blocknumber + noBlocks * (pheno - 1)) = corr(betaBLUPxLDGM{1},betaBLUPxUKBLD{1})^2;
        
        disp(T(blocknumber + noBlocks * (pheno - 1),:))
        
    end
    rows = noBlocks * (pheno - 1) + (1:noBlocks);
    phenoTable.phenotype(pheno) = pheno_name_array(pheno);
    phenoTable.pval_mse(pheno) = binocdf(sum(T.ldcov_mse(rows)>T.ldgm_mse(rows))-1,noBlocks,.5,'upper');
    phenoTable.pval_weights(pheno) = binocdf(sum(T.ldgm_r2(rows)>T.ldcov_r2(rows))-1,noBlocks,.5,'upper');
    phenoTable.mean_ldcov_mse(pheno) = mean(T.ldcov_mse(rows));
    phenoTable.mean_ldgm_mse(pheno) = mean(T.ldgm_mse(rows));
    phenoTable.mean_ldcov_weights_r2(pheno) = mean(T.ldcov_r2(rows));
    phenoTable.mean_ldgm_weights_r2(pheno) = mean(T.ldgm_r2(rows));
    
end
writetable(T,save_path);
writetable(phenoTable,phenoTable_save_path);


% ​
% ​    prediction_correlation_P(1) = qfnP(betaBLUPxLDGM_cat,betaBLUPxLDcov_cat,P_cat,whichSNPs_cat) / ...
%         sqrt(qfnP(betaBLUPxLDGM_cat,betaBLUPxLDGM_cat,P_cat,whichSNPs_cat) ...
%         * qfnP(betaBLUPxLDcov_cat,betaBLUPxLDcov_cat,P_cat,whichSNPs_cat));
%
%     prediction_correlation_P(2) = qfnP(betaBLUPxLDGM_cat,betaBLUPxUKBLD_cat,P_cat,whichSNPs_cat) / ...
%         sqrt(qfnP(betaBLUPxLDGM_cat,betaBLUPxLDGM_cat,P_cat,whichSNPs_cat) ...
%         * qfnP(betaBLUPxUKBLD_cat,betaBLUPxUKBLD_cat,P_cat,whichSNPs_cat));
%
%     prediction_correlation_P(3) = qfnP(betaBLUPxLDcov_cat,betaBLUPxUKBLD_cat,P_cat,whichSNPs_cat) / ...
%         sqrt(qfnP(betaBLUPxUKBLD_cat,betaBLUPxUKBLD_cat,P_cat,whichSNPs_cat) ...
%         * qfnP(betaBLUPxLDcov_cat,betaBLUPxLDcov_cat,P_cat,whichSNPs_cat));