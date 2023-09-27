function [mus_samples,Cs_samples] = SitePDFKnowledge_MusicX_fun(data, johnparas, jmtype, sofv)
% some parameters
depthdata = data(:,1); ydata = data(:,2:end);
[nrow, ncol] = size(ydata);
T = 15000; tb = 5000;
mu0 = zeros(ncol,1); C0 = diag((10^4)*ones(ncol,1),0);
degree_free = ncol+1; shape_para = 0.5; A = 100;
sof = sofv; % unit: m
% Normalization by Johnson distribution system
xdata = ydata;
for j = 1:ncol
    jparas = johnparas(j,:);
    for i = 1:nrow
        if ~isnan(ydata(i,j))
            y = ydata(i,j);
            xdata(i,j) = Johnson_member_yTOx(y, jmtype{j}, jparas);
        end
    end
end
xdata_col_ = xdata'; xdata_col = xdata_col_(:);
% correlation matrix
for i = 1:nrow
    hi = depthdata(i);
    for j = 1:nrow
        hj = depthdata(j);
        dh = abs(hi-hj);
        rho(i,j) = exp(-0.5*dh/sof);
    end
end
% Gibbs sampling
% step-1: initialization
it = 1;
rng(it)
newdata(:,:,it) = xdata;
for i = 1:ncol
    indx_nan = find(isnan(xdata(:,i)));
    newdata(indx_nan,i,it) = 0;
end
mus(:,it) = mu0 + sqrtm(C0)*randn(ncol,1);
for j = 1:ncol
    a_vector(j, it) = 1/gamrnd(0.5, A^2);
end
scale_matrix = 2*(degree_free-ncol+1)*diag(1./a_vector(:, it),0);
Cs(:,:,it) = iwishrnd(scale_matrix, degree_free);
inv_Cs(:,:,it) = inv(Cs(:,:,it));
vector_1 = ones(nrow,1);
x_col_ = newdata(:,:,it)'; x_col = x_col_(:);
x_mat = newdata(:,:,it)';
while it <= T-1
    % step-2: draw mu_s samples
    mu_part1 = inv(inv(C0)+(vector_1'*inv(rho)*vector_1)*inv_Cs(:,:,it));
    mu_part2 = kron((vector_1'*inv(rho)),inv_Cs(:,:,it));
    mus_mu = mu_part1*mu_part2 * x_col;
    mus_covarance = mu_part1;
    mus(:,it+1) = mus_mu + sqrtm(mus_covarance)*randn(ncol,1);
    % step-3: draw Cs samples
    Cs_part2 = (x_mat-mus(:,it+1)*vector_1')*inv(rho)*...
        (x_mat-mus(:,it+1)*vector_1')';
    Cs_scale_matrix = scale_matrix + Cs_part2;
    e = (1e-6)*ones(size(Cs_scale_matrix)); 
    % e: small scale perturbation. To ensure the positive matrix
    Cs_degree_free = nrow+ncol+1;
    Cs(:,:,it+1) = iwishrnd(Cs_scale_matrix+e, Cs_degree_free);
    inv_Cs(:,:,it+1) = inv(Cs(:,:,it+1));
    % step-4: draw a_vector samples
    shape_para = (ncol+2)/2;
    for j = 1:ncol
        scale_para = 1/A^2 + 2*inv_Cs(j,j,it+1);
        a_vector(j, it+1) = 1/gamrnd(shape_para , 1/scale_para);
    end
    scale_matrix = 2*(degree_free-ncol+1)*diag(1./a_vector(:, it+1));
    % step-5: draw x_u samples
    newdata(:,:,it+1) = newdata(:,:,it);
    vector_1_mus_kron = kron(vector_1, mus(:,it+1));
    rho_Cs_kron = kron(rho, Cs(:,:,it+1));
    indx_u = find ( isnan(xdata_col) ); % Find the locations of NaNs.
    indx_o = find ( ~isnan(xdata_col) );
    if indx_u ~= 0
        x_u = x_col(indx_u);
        x_o = x_col(indx_o);
        vector_1_mus_kron_u = vector_1_mus_kron(indx_u);
        vector_1_mus_kron_o = vector_1_mus_kron(indx_o);
        rho_Cs_kron_u = rho_Cs_kron(indx_u, indx_u);
        rho_Cs_kron_o = rho_Cs_kron(indx_o, indx_o);
        rho_Cs_kron_uo = rho_Cs_kron(indx_u, indx_o);
        rho_Cs_kron_ou = rho_Cs_kron_uo';
        % conditioning
        mu_xu = vector_1_mus_kron_u + rho_Cs_kron_uo*...
            inv(rho_Cs_kron_o)*(x_o-vector_1_mus_kron_o);
        Cs_xu = rho_Cs_kron_u - rho_Cs_kron_uo*...
            inv(rho_Cs_kron_o)*rho_Cs_kron_ou;
        % sampling xi_u
        x_u = mu_xu + sqrtm(Cs_xu)*randn(length(indx_u),1);
        % replacing
        x_col_ = newdata(:,:,it+1)'; x_col = x_col_(:);
        x_col(indx_u) = x_u;
        newdata(:,:,it+1) = (reshape(x_col, ncol, nrow))';
    end
    x_col_ = newdata(:,:,it+1)'; x_col = x_col_(:);
    x_mat = newdata(:,:,it+1)';
    it = it + 1;
    rng(it)
end
mus_samples = mus(:, tb+1:end);
Cs_samples = Cs(:, :, tb+1:end);
%inv_Cs_samples = inv_Cs(:, :, tb+1:end);
end