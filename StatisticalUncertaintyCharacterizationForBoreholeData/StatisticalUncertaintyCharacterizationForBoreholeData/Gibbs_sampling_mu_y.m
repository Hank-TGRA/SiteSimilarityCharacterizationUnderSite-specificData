function [mu_y]=Gibbs_sampling_mu_y(c,y,n,mu_mu0,sigma2_mu0,sigma2_y0)
% MCMC simulation regarding mean equivalent samples using Gibbs algorithm.
% --
% obtain the inverse of the matrix c.
inv_c = inv(c);
sub_eq1=0;
sub_eq2=0;
for i=1:1:n
    for j=1:1:n
        sub_eq1=inv_c(i,j)*(y(i)+y(j))+sub_eq1;
        sub_eq2=inv_c(i,j)+sub_eq2;
    end
end
% Calculate the distribution parameters according to the basic parameters.
mu_mu=(2*mu_mu0*sigma2_y0+sigma2_mu0*sub_eq1)/(2*(sigma2_y0+sigma2_mu0*sub_eq2));
sigma2_mu=(1/sigma2_mu0+sub_eq2/sigma2_y0)^(-1);
% Sampling the equivalent samples of mean according to the distribution parameters.
mu_y=random('Normal',mu_mu,sqrt(sigma2_mu));
end
