function y = kSystem(X,i,syspra,repN)

% X is the design matrix, X=[x1; x2; ... xm]
% x1, x2, .... xm, are m design points; each one is row vector
% i == 0, run all systems; otherwise, only run system i
% syspra, the parameters of systems
% repN, how many repeated observations

m = size(X,1); % number of design points
Z = [ones(m,1) X];

beta0 = syspra.beta0; % 1 x k matrix, k is system number
beta = syspra.beta;   % d x k matrix, d is dimension of design points
theta = [beta0; beta];% (1+d) x k

k = length(beta0);
homo = syspra.homo;

if homo == 1     % output (homo) noisy observations
    sigma2 = syspra.sigma2;
    sigma2 = repmat(sigma2, m,1);
elseif homo == 0 % output (heter) noisy observations
    sigma2 = 100 * (Z*theta).^2;
end

% randn, generate standard normal
if i == 0 % rum all systems
    y = repmat(Z*theta,[1,1,repN]) + randn(m,k,repN).*repmat(sqrt(sigma2),[1,1,repN]);
else % run system i
    y = repmat(Z*theta(:,i),[1,1,repN]) + randn(m,1,repN).*repmat(sqrt(sigma2(:,i)),[1,1,repN]);
end

end