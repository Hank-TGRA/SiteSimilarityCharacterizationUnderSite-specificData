function [sigma2_y]=Gibbs_sampling_sigma2_y(c,y,n,alpha_sigma0,beta_sigma0,mu_y)
% MCMC simulation regarding var equivalent samples using Gibbs algorithm.
% --
% Calculate the parameters of Inverse Gamma distribution
alpha_sigma=n/2+alpha_sigma0;
beta_sigma=beta_sigma0+0.5*((y-mu_y)'/c)*(y-mu_y);
% Samping var equivalent samples according to the distribution parameters.
sigma2_y=1/gamrnd(alpha_sigma,1/beta_sigma);
end