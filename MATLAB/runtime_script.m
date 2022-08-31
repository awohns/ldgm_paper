% This script runs + times equivalent precision matrix + correlation matrix 
% operations 
% 
% Use it to reproduce results of Figure 4a
% 
% Updated 8/25/22

clear;clc

jobnumber = str2num(getenv('SGE_TASK_ID'));
if isempty(jobnumber);jobnumber=1;end

addpath(genpath('/broad/oconnor/trees/ldgm/'));

ldgms_dir = '/broad/oconnor/trees/LDGMs_aug2022/';
save_file = '/broad/oconnor/trees/results_081522/runtime_082522.txt';

missingness = 0.1; %proportion of missing SNPs

reps = 1;
pop = 'AFR';
nn = 1e4;
h2 = 0.1;

% Load precision matrix
precision_fnames = dir([ldgms_dir,'1kg_chr*.',pop,'.edgelist']);
pm_file = [precision_fnames(jobnumber).folder, '/', precision_fnames(jobnumber).name];
P = readedgelist(pm_file);

% Correlation matrix
SNPs = diag(P) ~= 0;
P = P(SNPs,SNPs);
R = inv(full(P));

time = zeros(8,2);
for rep = 1:reps

    % Simulate sumstats
    
    if 1
        % 'Allele frequency'
        AF = 0.5 * full(diag(P) ~= 0);

        tic;
        [sumstats, ~, ~, beta] = simulateSumstats(...
            nn, {AF},...
            'precisionMatrices',{P},...
            'heritability',h2);
        time(3,2) = time(3,2) + toc;
        beta = beta{1};
        
        tic;
        simulateSumstats(...
            nn, {AF},...
            'correlationMatrices',{R},...
            'heritability',h2);
        time(3,1) = time(3,1) + toc;

    end

    % Missing data
    whichSNPs = (rand(length(P),1) < 1 - missingness);
    sumstats{1} = sumstats{1}(whichSNPs,:);
    Rs = R(whichSNPs,whichSNPs);
    noSNPs = sum(whichSNPs);

    % Z scores
    z = sumstats{1}.Z_deriv_allele;
    

    % Multiply by R
    if 1
        tic;
        x0 = Rs * z;
        time(1,1) = time(1,1) + toc;
        
        % Divide by P
        tic;
        x = precisionDivide(P, z, whichSNPs);
        time(1,2) = time(1,2) +toc;
        
    end
    
    % Divide by R
    if 1
        tic;
        xx = Rs \ z;
        time(2,1) = time(2,1) + toc;
        
        % Multiply by P
        tic;
        x = precisionMultiply(P, z, whichSNPs);
        time(2,2) = time(2,2) + toc;
        
    end
    
    % Likelihood given beta using R
    if 1
        
        tic;
        mu = Rs * randn(size(z));
        ll0 = mvnpdf(z, ...
            mu, 1/nn * Rs);
        
          % This approach is slower
%         mu = R(whichSNPs,:) * beta;
%         x = sumstats0 - mu;
%         ll0 = -1/2 * (log(det(R(whichSNPs,whichSNPs))) - log(nn)) - ...
%             1/2 * x' * (R(whichSNPs,whichSNPs) \ x);
        time(4,1) = time(4,1) + toc;
        
        
        % Likelihood given beta using P
        tic;
        mu = precisionDivide(P,randn(size(z)),whichSNPs);
        sigmasq = zeros(sum(whichSNPs),1);
        ll = GWASlikelihood(z - mu, sigmasq, ...
            P, nn, whichSNPs);
        time(4,2) = time(4,2) + toc;
        
    end
    
    % Likelihood given Sigma using R
    if 1
        Sigma = mean(beta.^2) * ones(sum(whichSNPs),1);
        tic;
        ll0 = log(mvnpdf(z, zeros(size(z)), ...
            1/nn * Rs + ...
            Rs * diag(Sigma) * Rs));
        time(5,1) = time(5,1) + toc;
        
        % Likelihood given Sigma using P
        tic;
        ll = GWASlikelihood(z, Sigma, P, nn, whichSNPs);
        time(5,2) = time(5,2) + toc;
        
    end
    
    % Compute the determinant
    if 1
        tic;
        A = chol(Rs);
        logDetR0 = 2*sum(log(diag(A)));
        
        % Slower
%         logDetR0 = log(det(R(whichSNPs,whichSNPs)));
        time(6,1) = time(6,1) + toc;
        
        tic;
        logDetR = logDeterminant(P,whichSNPs,zeros(size(z)));
        time(6,2) = time(6,2) + toc;
    end
    
    % BLUP
    if 1
        betaCov = mean(beta.^2);
        tic;
        betaBLUP0 = BLUPxldcov({Rs},{true(noSNPs,1)},sumstats,betaCov,nn);
        time(7,1) = time(7,1) + toc;

        tic;
        betaBLUP = BLUPxldgm({P},{whichSNPs},sumstats,betaCov,nn);
        time(7,2) = time(7,2) + toc;
        
    end
            
    % Likelihood gradient
    if 1
        tic;
        RSinv = inv(Rs.*Sigma' + 1/nn * eye(sum(whichSNPs)));
        RSinv_sumstats = RSinv * z/sqrt(nn);
        grad0 = -1/2 * sum(RSinv .* Rs)' + ...
            1/2 * sum(RSinv_sumstats .* RSinv_sumstats')';
        time(8,1) = time(8,1) + toc;
        
        tic;
        grad = GWASlikelihoodGradient(z/sqrt(nn),Sigma,P,nn,eye(length(Sigma)),whichSNPs,1);
        time(8,2) = time(8,2) + toc;
    end
    disp(time)
end
time = time/reps;

file = fopen(save_file,'a');
fprintf(file,'%f %f\n',time');
fclose(file);







