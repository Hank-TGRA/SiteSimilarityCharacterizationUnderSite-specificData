clc; clear; close;
datasetname = 'BoreholeData_LI_selected 31 sites';
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
    sitedata0 = data(indx,2:end);
    nonanindx = find(~isnan(sitedata0(:,2)));
    sitedata = sitedata0(nonanindx,:);
    k = k+1;
    mean_site(k,1) = mean(sitedata(:,2));
    var_site(k,1) = var(sitedata(:,2),1);
    [sof_site(k,1), ~,Lcn(k,1), Lindx(k,1),cnmin(k,1), cnminindx(k,1)] = MLE_SOF(sitedata, ACF);
    
end
initial_statistics = [site_id,mean_site,var_site,sof_site];
toc;
