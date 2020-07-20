function hstar = find_h_het_PCSmin(pra,hL,hU)
% This function calculates the constant h_min for heteroscedastic errors.
% With proper value of Xmin, h is the desired h_min in the paper.

% X is random vector, column vector
% X_i ~ Unif[0,1]

alpha = pra.alpha;
Xmin = pra.Xmin;

% Waypoints is used to avoid the numerical error in integral when degree of 
% degree of chi square distribution is too high!
% het, df = n0-1
if length(Xmin) > 1
    n0 = pra.n0; % initial sample size
    chix1 = fzero(@(chix) chi2pdf(chix,n0-1) - chi2pdf(n0-1,n0-1)*(1e-7) ,[0 n0-1]);
    chix2 = fzero(@(chix) chi2pdf(chix,n0-1) - chi2pdf(n0-1,n0-1)*(1e-7) ,[n0-1 (n0-1)*3]);
    Waypoints_het = [chix1 chix2];
    pra.Waypoints_het = round(Waypoints_het);
end

fh2 = @(h) fun_outer(Xmin',h,pra);
hstar = fzero(@(h) fh2(h) - (1-alpha),[hL hU]);   
end

%%% Given h, X, compute the outer function
function f_out = fun_outer(X,h,pra)

m = pra.m; % design points
n0 = pra.n0; % initial sample size

if length(X) > 1
    Waypoints_het = pra.Waypoints_het;
    f_out = integral(@(y) f_inner(y,X,h,pra) .* (m * chi2pdf(y,n0-1) .* (1-chi2cdf(y,n0-1)).^(m-1)), 0,Inf, 'Waypoints',Waypoints_het, 'AbsTol',1e-1,'RelTol', 1e-1);
else
    f_out = integral(@(y) f_inner(y,X,h,pra) .* (m * chi2pdf(y,n0-1) .* (1-chi2cdf(y,n0-1)).^(m-1)), 0,Inf);
end
end

%%% Given h, X, y, compute the inner function
function f_in = f_inner(yv,X,h,pra)

L = length(yv);
f_in = zeros(1,L);

k = pra.k; % system number
m = pra.m; % design points
Z = pra.Z;
n0 = pra.n0; % initial sample size

if length(X) > 1
    Waypoints_het = pra.Waypoints_het;
end

for l = 1:L
    y = yv(l);
    fun1 = @(x) ...
        normcdf(1./sqrt(1/y + 1./x) * h/sqrt((n0-1)*([1 X']/(Z'*Z)*[1; X])))...
        .* (m * chi2pdf(x,n0-1) .* (1-chi2cdf(x,n0-1)).^(m-1));
    if length(X) > 1
        f_in(l) = integral(fun1,0,Inf, 'Waypoints',Waypoints_het, 'AbsTol',1e-1,'RelTol', 1e-1);
    else
        f_in(l) = integral(fun1,0,Inf);
    end
    f_in(l) = f_in(l)^(k-1);
end
end