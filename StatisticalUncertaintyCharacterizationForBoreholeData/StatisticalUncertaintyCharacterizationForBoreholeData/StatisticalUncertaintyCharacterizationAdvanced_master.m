clc; clear; close all;
% Loading Borehole data
% Note: data is a n*3 matrix, where n represents the number of data
% samples, 1st column contains the site id info, 2nd contains the depth
% data, 3rd contains the geo-material parametric data.
BoreholeDataName = 'BoreholeData_LI_Onsoy site';
InitialStatisticName = 'InitialStatistics_LI_Onsoy site';
filesavepath = 'D:\CodeTest\SUCA_LI_Onsoy';
datatype = 1; % 0: original data; 1: ln transformed data;
ACF = 1; % Select the auto-correlation function
BaiscParameters=[1.04, 0.408; 
                 2.252, 0.217;
                 2.031, 8.016];
% In BaiscParameters, 1st row contains the basic parameters required by
% posterior mean distribution, namely the mean of mean, and var of mean;
% 2nd contains the basic parameters required by the posterior var,
% calculated by the specific euqation in Han et al. (2022); 3rd row
% contains the basic parameters required by posterior sof, namely the mean
% of sof, and the var of sof.
%
Borehole_data = load([BoreholeDataName, '.txt']);
Initial_statistic = load([InitialStatisticName, '.txt']);
filepath_ACF = [filesavepath, '\ACF', num2str(ACF), '_'];
siteid = unique(Borehole_data(:,1)); numsite = length(siteid);
for i = 1:numsite
    id = siteid(i);
    disp(['Analysis site ID: ', num2str(id), ' (ACF-', num2str(ACF), ')']);
    indx = find(Borehole_data(:,1)==id);
    sitedata0 = Borehole_data(indx,2:end);
    nonanindx = find(~isnan(sitedata0(:,2)));
    if length(nonanindx)>=6
        sitedata = sitedata0(nonanindx,:);
        indx = find(Initial_statistic(:,1)==id);
        initialsta = Initial_statistic(indx,2:end);
        % Initial values for var and sof: [initial_var, initial_sof]
        Initial_var_and_theta = [initialsta(1), initialsta(2)];
        % Basic parameter vector of posterior distribution regarding mean:
        % [mean_of_mean, var_of_mean]
        PostPDF_Para_mu = BaiscParameters(1,:);
        % Basic parameter vector of posterior distribution (inverse Gamma 
        % function) regarding var: [alpha_sigma0, beta_sigma0]
        PostPDF_Para_var = BaiscParameters(2,:);
        % Basic parameter vector of posterior distribution regarding sof:
        % [mean_of_sof, var_of_sof]
        PostPDF_Para_theta = BaiscParameters(3,:);
        % -----------------------------------------------------;
        % call the function for statistical uncertainty characterization
        [equivlent_samples] =  StatisticalUncertaintyCharacterization(sitedata, Initial_var_and_theta,...
            PostPDF_Para_mu, PostPDF_Para_var, PostPDF_Para_theta, ACF);
        % Output the simulated equivalent samples in .txt format, number of 
        % data samples is 10000.
        % Note: Program generates 12000 equivalent samples, while the front
        % 2000 equivalent samples are regarded as the "burn-in period" 
        % samples and removed.
        filenameoutput1 = ['equivlent_samples-site-', num2str(id), '.txt'];
        fid_samples=fopen([filepath_ACF, filenameoutput1],'wt');
        fprintf(fid_samples,'%s %s %s\n','Mean','Variance','AutoCorreDistance');
        fprintf(fid_samples,'%f %f %f\n',equivlent_samples');
        fclose(fid_samples);
        % obtian the quantiles
        sample_2_5 = prctile(equivlent_samples,2.5,1);
        sample_5 = prctile(equivlent_samples,5,1);
        sample_10 = prctile(equivlent_samples,10,1);
        sample_50 = prctile(equivlent_samples,50,1);
        sample_90 = prctile(equivlent_samples,90,1);
        sample_95 = prctile(equivlent_samples,95,1);
        sample_97_5 = prctile(equivlent_samples,97.5,1);
        sample_mean = mean(equivlent_samples);
        sample_CI95 = sample_97_5-sample_2_5;
        sample_CI90 = sample_95-sample_5;
        sample_CI80 = sample_90-sample_10;
        sample_mu_statistics(i,:) = [id, sample_5(1), sample_10(1), sample_50(1),...
            sample_90(1), sample_95(1), sample_mean(1), sample_CI95(1), sample_CI90(1), ...
            sample_CI80(1)];
        sample_var_statistics(i,:) = [id, sample_5(2), sample_10(2), sample_50(2),...
            sample_90(2), sample_95(2), sample_mean(2), sample_CI95(2), sample_CI90(2), ...
            sample_CI80(2)];
        sample_sof_statistics(i,:) = [id, sample_5(3), sample_10(3), sample_50(3),...
            sample_90(3), sample_95(3), sample_mean(3), sample_CI95(3), sample_CI90(3), ...
            sample_CI80(3)];
        %
        if datatype == 1
            mu_lower = exp(sample_10(1)+sample_10(2)/2);
            mu_upper = exp(sample_90(1)+sample_90(2)/2);
            mu_updated = exp(sample_mean(1)+sample_mean(2)/2);
            sd_updated = sample_mean(1)*sqrt(exp(sample_mean(2)-1));
            sof_updated = sample_mean(3);
            sample_rf_params(i,:) = [id, mu_updated, sd_updated, sof_updated, mu_lower, mu_upper];
        elseif datatype == 0
            mu_lower = sample_10(1);
            mu_upper = sample_90(1);
            mu_updated = sample_mean(1);
            sd_updated = sqrt(sample_mean(2));
            sof_updated = sample_mean(3);
            sample_rf_params(i,:) = [id, mu_updated, sd_updated, sof_updated, mu_lower, mu_upper];
        end
    else
        sample_mu_statistics(i,:) = [id, NaN*ones(1,9)];
        sample_var_statistics(i,:) = [id, NaN*ones(1,9)];
        sample_sof_statistics(i,:) = [id, NaN*ones(1,9)];
        sample_rf_params(i,:) = [id, NaN*ones(1,5)];
        continue;
    end
end
%
filenameoutput2 = 'equivlent_mu_statistics.txt';
fid2=fopen([filepath_ACF, filenameoutput2],'wt');
fprintf(fid2,'site_id, p5, p10, p50, p90, p95, mean, CI95, CI90, CI80\n');
fprintf(fid2,'%g, %f, %f, %f, %f, %f, %f, %f, %f, %f\n',sample_mu_statistics');
fclose(fid2);
%
filenameoutput3 = 'equivlent_var_statistics.txt';
fid2=fopen([filepath_ACF, filenameoutput3],'wt');
fprintf(fid2,'site_id, p5, p10, p50, p90, p95, mean, CI95, CI90, CI80\n');
fprintf(fid2,'%g, %f, %f, %f, %f, %f, %f, %f, %f, %f\n',sample_var_statistics');
fclose(fid2);
%
filenameoutput4 = 'equivlent_sof_statistics.txt';
fid2=fopen([filepath_ACF, filenameoutput4],'wt');
fprintf(fid2,'site_id, p5, p10, p50, p90, p95, mean, CI95, CI90, CI80\n');
fprintf(fid2,'%g, %f, %f, %f, %f, %f, %f, %f, %f, %f\n',sample_sof_statistics');
fclose(fid2);
%
filenameoutput5 = 'CRF_Params.txt';
fid5=fopen([filepath_ACF, filenameoutput5],'wt');
fprintf(fid5,'site_id, mu_updated, sd_updated, sof_updated, mu_lower, mu_upper\n');
fprintf(fid5,'%g %f %f %f %f %f\n',sample_rf_params');
fclose(fid5);