clc; clear; close;
datasetname = 'BackgroundData_LI_Clay_10_7490';
ACF = 1;
% ---------------------------------------
tic;
data = load([datasetname, '.txt']);
[n,~] = size(data);
site_id = unique(data(:,1));
numsite = length(site_id);
k = 0;
for i = 1:numsite
    indx = find(data(:,1) == site_id(i));
    sitedata0 = data(indx,[2,3]);
    nonanindx = find(~isnan(sitedata0(:,2)));
    sitedata = sitedata0(nonanindx,:);
    numsamples = length(nonanindx);
    if numsamples >= 6
        k= k+1;
        site_id_use(k,:) = [k, site_id(i), numsamples];
        mean_data(k,1) = mean(sitedata(:,2));
        var_data(k,1) = var(sitedata(:,2),1);
        [sof_data(k,1), ~,Lcn(k,1), Lindx(k,1),cnmin(k,1), cnminindx(k,1)] = MLE_SOF(sitedata, ACF);
    end
end
%
indx = find(~isnan(sof_data(:)));
sof_data_NaN_remove = sof_data(indx);
%
mean_mean = mean(mean_data);
mean_media = prctile(mean_data,50);
mean_var = var(mean_data, 1);
var_mean = mean(var_data);
var_media = prctile(var_data,50);
var_var = var(var_data, 1);
sof_mean = mean(sof_data_NaN_remove);
sof_media = prctile(sof_data_NaN_remove,50);
sof_var = var(sof_data_NaN_remove,1);
%
alpha_sigma0 = var_mean^2/var_var+2;
beta_sigma0 = var_mean*(alpha_sigma0-1);
%
fprintf('Mean_of_Mean, mean_mean = %.3f; Var_of_Mean, mean_var = %.3f; mean_media = %.3f\n',mean_mean, mean_var,mean_media);
fprintf('Mean_of_Var, var_mean = %.3f; Var_of_Var, var_var = %.3f; var_media = %.3f\n',var_mean, var_var,var_media);
fprintf('alpha_sigma0 = %.3f; beta_sigma0 = %.3f\n',alpha_sigma0, beta_sigma0);
fprintf('Mean_of_SOF, sof_mean = %.3f; Var_of_SOF, sof_var = %.3f; sof_media = %.3f\n',sof_mean, sof_var,sof_media);
toc;
