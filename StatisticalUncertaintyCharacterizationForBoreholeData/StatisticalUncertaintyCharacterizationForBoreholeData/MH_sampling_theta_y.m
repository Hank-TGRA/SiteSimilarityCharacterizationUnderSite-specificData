function sof=MH_sampling_theta_y(c,y,n,sigma2_y,mu_y,sigma2_theta0,mu_theta0,theta_y0)
% assume that sof follows lognormal distribution
var_ln_sof = log(1+sigma2_theta0/(mu_theta0^2)); 
mu_ln_sof = log(mu_theta0)-0.5*var_ln_sof;
ln_sof_0 = log(theta_y0);
%
mu_y_vec = mu_y*ones(size(y));
% MCMC simulation regarding SOF equivalent samples using MH algorithm.
% postprior distribution of theta_y
syms x;
sub_eq1=1/sqrt(((2*pi*sigma2_y)^n)*det(c));
sub_eq2=exp(-(1/(2*sigma2_y))*((y-mu_y_vec)'/c)*(y-mu_y_vec));
% sub_eq3=1/sqrt(2*pi*sigma2_theta0);
% sub_eq4=exp(-((x-mu_theta0)^2)/(2*sigma2_theta0));
sub_eq3=1/(sqrt(2*pi*var_ln_sof)*theta_y0);
sub_eq4=exp(-((x-mu_ln_sof)^2)/(2*var_ln_sof));
f=sub_eq1*sub_eq2*sub_eq3*sub_eq4;
% M-H MCMC for sampling theta_y, proposal pdf is norm distribution
proposal_COV_theta = sqrt(var_ln_sof)/mu_ln_sof;
ln_sof_1=random('Normal',ln_sof_0,abs(ln_sof_0*proposal_COV_theta));
% "Accept/reject" procedure for samples with values > 0
% if ln_sof_1>0
% Calculate the acceptance rate ra
pdf_ln_sof_1=eval(subs(f,x,ln_sof_1));
pdf_ln_sof_0=eval(subs(f,x,ln_sof_0));
pdf_pro_ln_sof_0 = normpdf(ln_sof_0,ln_sof_1,abs(ln_sof_1*proposal_COV_theta));
pdf_pro_ln_sof_1 = normpdf(ln_sof_1,ln_sof_0,abs(ln_sof_0*proposal_COV_theta));
ra=(pdf_ln_sof_1*pdf_pro_ln_sof_0)/(pdf_ln_sof_0*pdf_pro_ln_sof_1);
% Generate the random number within the range of [0,1]
u=rand(1);
% Determine the value of a new sample according to the magnitude
% relationship between ra and u.
if ra>u
    ln_sof = ln_sof_1;
else
    ln_sof = ln_sof_0;
end
% else
%     theta_y=theta_y0;
% end
sof = exp(ln_sof);
end