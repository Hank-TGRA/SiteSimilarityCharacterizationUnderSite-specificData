function [final_theta,L_theta_vs_theta, Lcn, Lindx, cnmin, cnminindx] = MLE_SOF(borehole_data,ACF)
% This program is used to estimate the vertical SOF of a geo-material
% parameter obtained by the field investigation using the maximum
% likelihood method (MLM).
% --
% Note:
% According to relevant researches, the geo-material parameters can be
% assumed to follow the normal or lognormal distribution models. Since some
% geo-material parameters are non-negative, the natural lognorithm
% transform is applied to these non-negative parameters.
% --
% Meaning of ACF code:
% ACF=1: Markovian (single exponential) model
% ACF=2: Spherical model
% ACF=3: Second-order Markov model
% ACF=4: Gaussian (squared exponential) model
% --
% Loading the borehole data
data = borehole_data;
% Basic characteristics of data
[n, ~]= size(data);
depth = data(:,1);
y = data(:,2);
% Mean and standard deviation of geo-material parameters
mu_data = mean(y);
var_data = var(y,1);
% Set the number of iterative calculation
iteration = 200000;
% Define the vector for storaging data
L_theta = zeros(iteration,1);
theta = zeros(iteration,1);
% Set the initial sof value
thera0 = 0.001;
% Set the step of sof for iterative calculation
d_theta = 0.0001;
% Perform the iterative calculation
cn = zeros(iteration,1);
for i = 1:iteration
    theta(i) = thera0 + (i-1)*d_theta;
    corr_matrix  = Correlation_matrix(depth, n, theta(i), ACF);
    covar_matrix = corr_matrix.*var_data;
    cn(i) = rcond(covar_matrix);
    ii = 0;
    if cn(i) < 0.1
        L_theta(i) = NaN;
    else
        ii = ii+1;
        indxc(ii) = i;
        sub_equ = y-mu_data;
        L_theta(i) = -(0.5*n)*log(2*pi())-0.5*log(det(covar_matrix))...
            -0.5*((sub_equ'/covar_matrix)*sub_equ);
    end
end
% Visualization of results: plot the curve of SOF values versus likelihood function values
figure(1)
if ii ~= 0
    plot(theta(indxc), L_theta(indxc));
else
    plot(theta, L_theta);
end
xlabel('SOF value');  ylabel('likelihood function value');  
hold off
% Figure out the SOF value correponding to the maximum likelihood function value.
[~,Lindx]=max(L_theta);
if isnan(L_theta(Lindx))
    final_theta = NaN;
else
    final_theta = theta(Lindx);
end
Lcn = cn(Lindx);
% Print the most possible SOF value in the command window.
% fprintf('final_theta = %g\n', final_theta);
% Merge the SOF vector and the maximum likelihood function value vector to
% a matrix.
L_theta_vs_theta = [theta, L_theta];
% Identify the singular matrix
[cnmin, cnminindx] = min(cn);
end
