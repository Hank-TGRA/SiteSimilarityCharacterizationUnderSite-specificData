function [after_burn_sample] = StatisticalUncertaintyCharacterization(Input_data, Initial_var_and_theta,...
    PostPDF_Para_mu, PostPDF_Para_var, PostPDF_Para_theta, ACF)
% This program is used to characterize the statistical uncertainty of borehole data in geotechnical
% engineering field investigation.
% This program involves 3 statistical characteristics, including mean, var, SOF (scale of fluctuation).
% For the non-negative geo-material parameters, the natural logarithmic transformation should be 
% applied to them. This is because that the geo-material parameters are assumed to follow the normal
% and lognormal distribution models.
% Code meaning of ACF:
% ACF=1: Markovian (single exponential) model
% ACF=2: Spherical model
% ACF=3: Second-order Markov model
% ACF=4: Gaussian (squared exponential) model
% --
tic
% Progress bar
h=waitbar(0,'Calculating, please wait£¡');
% Loading data [depth data, parametric data]
data = Input_data;
y = data(:,2);
depth = data(:,1);
[n_instances, ~] = size(data);
% Set the initial values
sigma2_y0 =  Initial_var_and_theta(1);
theta_y0 = Initial_var_and_theta(2);
c = Correlation_matrix(depth,n_instances,theta_y0, ACF);
% a series of parameters of pdf
% for postprior distribution of mu_y
mu_mu0 = PostPDF_Para_mu(1);
sigma2_mu0 = PostPDF_Para_mu(2);
% for postprior distribution of sigma2_y
alpha_sigma0 = PostPDF_Para_var(1);
beta_sigma0 = PostPDF_Para_var(2);
% for postprior distribution of theta_y
mu_theta0 = PostPDF_Para_theta(1);
sigma2_theta0 = PostPDF_Para_theta(2);
% sampling procedure
ns=12000;  % number of sampling
sample=zeros(ns,3);
for i=1:1:ns
    rng(i, 'v4');
    % Set the progress bar
    str=['Task completion rate£º',num2str(i/ns*100),'%'];
    waitbar(i/ns,h,str);
    % sampling mu_y
    mu_y=Gibbs_sampling_mu_y(c,y,n_instances,mu_mu0,sigma2_mu0,sigma2_y0);
    % sampling sigma2_y
    sigma2_y0=Gibbs_sampling_sigma2_y(c,y,n_instances,alpha_sigma0,beta_sigma0,mu_y);
    % sampling theta_y
    theta_y0=MH_sampling_theta_y(c,y,n_instances,sigma2_y0,mu_y,sigma2_theta0,mu_theta0,theta_y0);
    % storing the samples
    sample(i,:)=[mu_y, sigma2_y0, theta_y0];
    % updating the correlation matrix using new theta_y
    c = Correlation_matrix(depth,n_instances,theta_y0, ACF);
end
% In 12000 samples, the front 2000 samples are seen as the 'burn-in period' samples. 
burn_period = (1/6)*ns; 
after_burn_sample = sample(burn_period+1:ns, :);
% Plot the histogram of mean equivalent samples
figure (1)
histogram(after_burn_sample(:,1));
legend('¦Ì_y')
xlabel('¦Ì_y');
ylabel('Number');
hold off
% Plot the histogram of var equivalent samples
figure (2)
histogram(after_burn_sample(:,2));
legend('¦Ò^2_y')
xlabel('¦Ò^2_y');
ylabel('Number');
hold off
% Plot the histogram of sof equivalent samples
figure (3)
histogram(after_burn_sample(:,3));
legend('¦È_y')
xlabel('¦È_y');
ylabel('Number');
hold off
close(h);
toc
end
