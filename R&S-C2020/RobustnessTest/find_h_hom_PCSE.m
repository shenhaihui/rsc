function hstar = find_h_hom_PCSE(pra,hL,hU,~)
% This function calculates the constant h for homoscedastic errors
% Using numerical integration

% X is random vector, column vector
% by default,      X_i ~ Unif[0,1]
% when specified,  X1 ~ N(0,1), X2 ~ N(0,1), X3 ~ N(0,1)
%                    truncated on [0,1], covariance = 0.5 for each pair

alpha = pra.alpha;
if nargin == 3
    pra.dist_type = 1;
else
    pra.dist_type = 2;
end
    
% Waypoints is used to avoid the numerical error in integral when degree of 
% chi square distribution is too high!
% hom, df = n0*m-1-d
d = pra.d; % dimension
if d > 1
    m = pra.m; % design points
    n0 = pra.n0; % initial sample size
    chix1 = fzero(@(chix) chi2pdf(chix,n0*m-1-d) - chi2pdf(n0*m-1-d,n0*m-1-d)*(1e-7) ,[0 n0*m-1-d]);
    chix2 = fzero(@(chix) chi2pdf(chix,n0*m-1-d) - chi2pdf(n0*m-1-d,n0*m-1-d)*(1e-7) ,[n0*m-1-d (n0*m-1-d)*2]);
    Waypoints_hom = [chix1 chix2];
    pra.Waypoints_hom = round(Waypoints_hom);
end

fh = @(h) f_h(h,pra);
if d == 1
    hstar = fzero(@(h) fh(h) - (1-alpha),[hL hU]); 
else
    options = optimset('TolX',0.001);
    hstar = fzero(@(h) fh(h) - (1-alpha+0.001),[hL hU],options);
end
end

%%% Adopt trapz function instead of integral function to evaluate the
%%% expectation of X, which is coarser but faster, when d > 1
function f3 = f_h(h,pra)

d = pra.d; % dimension

if d == 1
    f3 = integral(@(X) fun_outer1(X,h,pra),0,1);
elseif d ==2
    f1 = @(X1,X2) fun_outer2(X1,X2,h,pra);
    f2 = @(X1,X2) arrayfun(f1,X1,X2);
    f3 = integral2(f2,0,1,0,1); 
elseif d == 3
    if pra.dist_type == 1
        f1 = @(X1,X2,X3) fun_outer3(X1,X2,X3,h,pra);
        f2 = @(X1,X2,X3) arrayfun(f1,X1,X2,X3);
    else
        mu = [0.5;0.5;0.5];
        Sigma = [1 0.5 0.5;0.5 1 0.5; 0.5 0.5 1];
        f1 = @(X1,X2,X3) fun_outer3(X1,X2,X3,h,pra) * mvnpdf([X1;X2;X3],mu,Sigma)/0.0750;
        f2 = @(X1,X2,X3) arrayfun(f1,X1,X2,X3);
    end    
    X1g = 0:0.1:1;
    X2g = 0:0.1:1;
    X3g = 0:0.1:1;
    [X1m,X2m,X3m] = ndgrid(X1g,X2g,X3g);
    F = f2(X1m,X2m,X3m);
    f3 = trapz(X3g,trapz(X2g,trapz(X1g,F,1),2),3);   
elseif d == 5
    f1 = @(X1,X2,X3,X4,X5) fun_outer5(X1,X2,X3,X4,X5,h,pra);
    f2 = @(X1,X2,X3,X4,X5) arrayfun(f1,X1,X2,X3,X4,X5);
    X1g = 0:0.2:1;
    X2g = 0:0.2:1;
    X3g = 0:0.2:1;
    X4g = 0:0.2:1;
    X5g = 0:0.2:1;
    [X1m,X2m,X3m,X4m,X5m] = ndgrid(X1g,X2g,X3g,X4g,X5g); %7776
    F = f2(X1m,X2m,X3m,X4m,X5m);
    f3 = trapz(X5g,trapz(X4g,trapz(X3g,trapz(X2g,trapz(X1g,F,1),2),3),4),5);    
end
end


%%% Given h, X, compute the outer function
% d == 1
function f_out = fun_outer1(Xv,h,pra)
L = length(Xv);
f_out = zeros(1,L);
d = pra.d; % dimension
m = pra.m; % design points
n0 = pra.n0; % initial sample size
for l=1:L
    X = Xv(l);
    f_out(l) = integral(@(y) f_inner(y,X,h,pra) .* chi2pdf(y,n0*m-1-d),0,Inf);
end
end
% d == 2
function f_out = fun_outer2(X1,X2,h,pra)
d = pra.d; % dimension
m = pra.m; % design points
n0 = pra.n0; % initial sample size
Waypoints_hom = pra.Waypoints_hom;
X = [X1; X2];
f_out = integral(@(y) f_inner(y,X,h,pra) .* chi2pdf(y,n0*m-1-d),0,Inf, 'Waypoints',Waypoints_hom, 'AbsTol',1e-1,'RelTol', 1e-1);
end
% d == 3
function f_out = fun_outer3(X1,X2,X3,h,pra)
d = pra.d; % dimension
m = pra.m; % design points
n0 = pra.n0; % initial sample size
X = [X1; X2; X3];
Waypoints_hom = pra.Waypoints_hom;
f_out = integral(@(y) f_inner(y,X,h,pra) .* chi2pdf(y,n0*m-1-d),0,Inf, 'Waypoints',Waypoints_hom, 'AbsTol',1e-1,'RelTol', 1e-1);
end
% d == 5
function f_out = fun_outer5(X1,X2,X3,X4,X5,h,pra)
d = pra.d; % dimension
m = pra.m; % design points
n0 = pra.n0; % initial sample size
X = [X1; X2; X3; X4; X5];
Waypoints_hom = pra.Waypoints_hom;
f_out = integral(@(y) f_inner(y,X,h,pra) .* chi2pdf(y,n0*m-1-d),0,Inf, 'Waypoints',Waypoints_hom, 'AbsTol',1e-1,'RelTol', 1e-1);
end


%%% Given h, X, y, compute the inner function
function f_in = f_inner(yv,X,h,pra)

L = length(yv);
f_in = zeros(1,L);

k = pra.k; % system number
d = pra.d; % dimension
m = pra.m; % design points
Z = pra.Z;
n0 = pra.n0; % initial sample size
if d > 1
    Waypoints_hom = pra.Waypoints_hom;    
end

for l = 1:L
    y = yv(l);
    fun1 = @(x) ...
        normcdf(1./sqrt(1/y + 1./x) * h/sqrt((n0*m-1-d)*([1 X']/(Z'*Z)*[1; X])))...
        .* chi2pdf(x,n0*m-1-d);
    if d > 1
        f_in(l) = integral(fun1,0,Inf, 'Waypoints',Waypoints_hom, 'AbsTol',1e-1,'RelTol', 1e-1);
    else
        f_in(l) = integral(fun1,0,Inf);
    end
    f_in(l) = f_in(l)^(k-1);
end
end