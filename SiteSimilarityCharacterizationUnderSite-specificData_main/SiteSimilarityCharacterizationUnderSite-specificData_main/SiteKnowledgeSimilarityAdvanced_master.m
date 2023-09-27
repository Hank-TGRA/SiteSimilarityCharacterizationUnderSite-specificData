clc; clear; close;
filename_target_site = 'target_site_data';
selected_cp_sites_data = load('selected_cp_sites_data.txt');
selected_cp_sites_sof = load('selected_cp_sites_sof.txt');
johnparams = load('johnson_params_complete.txt');
jmtype = {'SU', 'SU','SU','SB','SB','SU','SU','SU','SU','SU'};
sofv_A = 1.809; % unit: m
filesavepath = 'D:\CodeTest\site_similarity_analysis';
% ----------------------------------------------------------
d = length(jmtype);
DataSiteA = load([filename_target_site, '.txt']);
siteid = unique(selected_cp_sites_data(:,1)); numsites = length(siteid);
smv_statistics_all = zeros(numsites,4);
for j = 1:numsites
    id = siteid(j);
    indx = find(selected_cp_sites_data(:,1)==id);
    DataSiteB = selected_cp_sites_data(indx,2:end);
    indx = find(selected_cp_sites_sof(:,1)==id);
    sofv_B = selected_cp_sites_sof(indx,2);
    if length(DataSiteA(1,2:end))~=d || length(DataSiteB(1,2:end))~=d
        error('Target site data, selected comparative site data, and Johnson params should have the same size! ')
    end
    %
    for i = 1:2
        if i == 1
            [mus_samplesA,Cs_samplesA]=SitePDFKnowledge_MusicX_fun(DataSiteA,johnparams,jmtype,sofv_A);
        elseif i == 2
            [mus_samplesB,Cs_samplesB]=SitePDFKnowledge_MusicX_fun(DataSiteB,johnparams,jmtype,sofv_B);
        end
    end
    [~, numsam] = size(mus_samplesA);
    smv = zeros(numsam,3);
    for i = 1:numsam
        [rho_matrixA, sds_vectorA]= ConvertCOVmatrixIntoRmatrixAndVarvector(Cs_samplesA(:,:,i));
        [rho_matrixB, sds_vectorB]= ConvertCOVmatrixIntoRmatrixAndVarvector(Cs_samplesB(:,:,i));
        %
        min_mu = min(min(mus_samplesA(:,i)), min(mus_samplesB(:,i)));
        mus_samplesA_re = mus_samplesA(:,i)-min_mu*ones(d,1); % re: revised
        mus_samplesB_re = mus_samplesB(:,i)-min_mu*ones(d,1);
        dynamicrange_mu = max(max(mus_samplesA_re), max(mus_samplesB_re));
        %
        min_sd = min(min(sds_vectorA), min(sds_vectorB));
        sds_samplesA_re = sds_vectorA - min_sd*ones(d,1);
        sds_samplesB_re = sds_vectorB - min_sd*ones(d,1);
        dynamicrange_sd = max(max(sds_samplesA_re), max(sds_samplesB_re));
        %
        rho_matrixA_re = rho_matrixA+ones(d,d);
        rho_matrixB_re = rho_matrixB+ones(d,d);
        dynamicrange_r = 2;
        %
        [smv(i,1), smp_mu(:,:,i)] = multissim(mus_samplesA_re, mus_samplesB_re,...
            'NumScales', 1, 'Sigma', 1.5,'DynamicRange', dynamicrange_mu);
        [smv(i,2), smp_sd(:,:,i)] = multissim(sds_samplesA_re, sds_samplesB_re,...
            'NumScales', 1, 'Sigma', 1.5,'DynamicRange', dynamicrange_sd); % cov vector
        [smv(i,3), smp_r(:,:,i)] = multissim(rho_matrixA_re, rho_matrixB_re,...
            'Sigma', 1.5,'DynamicRange', dynamicrange_r); % relationship matrix
    end
    smv_p50 = prctile(smv, 50); % median
    smv_p5 = prctile(smv, 5);
    smv_p95 = prctile(smv, 95);
    smv_mean = mean(smv,1);
    %
    smv_mu_mean = smv_mean(1);
    smv_sd_mean = smv_mean(2);
    smv_r_mean = smv_mean(3);
    %
    smv_mu_media = smv_p50(1);
    smv_sd_media = smv_p50(2);
    smv_r_media = smv_p50(3);
    %
    smv_statistics_all(j,:) = [id, smv_mu_mean, smv_sd_mean, smv_r_mean];
    % save
    fid1 = fopen([filesavepath,'\smv-', num2str(id), '.txt'], 'wt');
    fprintf(fid1, '%f %f %f\n', smv');
    fclose(fid1);
    fid2 = fopen([filesavepath,'\smv-', num2str(id), '-statistics.txt'], 'wt');
    fprintf(fid2, 'smv_mu_mean, smv_sd_mean, smv_r_mean\n');
    fprintf(fid2, '%f, %f, %f\n', smv_mu_mean, smv_sd_mean, smv_r_mean);
    fprintf(fid2, 'smv_mu_media, smv_sd_media, smv_r_media\n');
    fprintf(fid2, '%f, %f, %f\n', smv_mu_media, smv_sd_media, smv_r_media);
    fclose(fid2);
end
fid3 = fopen([filesavepath,'\smv-statistics-all.txt'], 'wt');
fprintf(fid3, 'site_id, smv_mu_mean, smv_sd_mean, smv_r_mean\n');
fprintf(fid2, '%g %f, %f, %f\n', smv_statistics_all');
fclose(fid2);
%% Covert the covariance matrix into correlation coefficient matrix and var vector
function [rmatrix, sdvector] = ConvertCOVmatrixIntoRmatrixAndVarvector(covmatrix)
varvector = diag(covmatrix);
sdvector = sqrt(varvector);
[~, d] = size(covmatrix);
for i=1:d
    for j=1:d
        sigma12 = sdvector(i)*sdvector(j);
        rmatrix(i,j) = covmatrix(i,j)/sigma12;
    end
end
end