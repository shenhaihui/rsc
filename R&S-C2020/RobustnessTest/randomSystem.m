function y = randomSystem(X,i,trueoutput,repN, f_true)

% X is the design matrix, X=[x1; x2; ... xm]
% x1, x2, .... xm, are m design points; each one is row vector
% i == 0, run all systems; otherwise, only run system i
% trueoutput == 0, noisy observations; otherwise, true values
% repN, how many repeated observations

m = size(X,1); % number of design points
k = 5; % number of alternatives

%%% home case, always compute for the whole design matrix

if trueoutput == 0 % noisy observations (Hom + Equal variance)
    sigma = 1;
    if i == 0 % rum all systems
        y = repmat(f_true,[1,1,repN]) + randn(m,k,repN) * sigma;
    else % run system i
        y = repmat(f_true(:,i),[1,1,repN]) + randn(m,1,repN) * sigma;
    end
end

end