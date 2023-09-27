function [cm] = Correlation_matrix(y,n,theta_y,ACF)
% Calculate the spatial auto-correlation coefficient matrix
cm=zeros(n,n);
for i=1:1:n
    r1=y(i,1);
    for j=1:1:n
        r2=y(j,1);
        d = abs(r1-r2);
        switch ACF
            case 1
                % Markovian (single exponential) model
                cm(i,j)=exp(-2*d/theta_y);
            case 2
                % Spherical model
                if d <= theta_y
                    cm(i,j) = (1-d/theta_y);
                else
                    cm(i,j) = 0;
                end
            case 3
                % Second-order Markov model
                cm(i,j) = exp(-4*d/theta_y)*(1+4*d/theta_y);
            case 4
                % Gaussian (squared exponential) model
                cm(i,j) = exp(-pi()*(d/theta_y)^2);
        end
    end
    
end
end