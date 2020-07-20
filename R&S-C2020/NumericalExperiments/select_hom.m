function [thetahat, totalN] = select_hom(h,pra,syspra)

% system number
k = pra.k;
% dimension
d = pra.d;
% design points
m = pra.m;
Xd = pra.Xd;
Z = pra.Z;
% initial sample size
n0 = pra.n0;
% indifferent-zone parameter
delta = pra.delta;


%%%%%% Step 0 - solve h

%%%%%% Step 1
% run n0 batches for k systems
Y = kSystem(Xd,0,syspra,n0); % n0 batches

Ybar = mean(Y,3);
S2 = zeros(1,k);
for i=1:k
    theta = (Z'*Z)\Z'*Ybar(:,i);    
    SV = 0;
    for l = 1:n0
        SV = SV + sum((Y(:,i,l) - Z*theta).^2);
    end
    S2(i) = SV/(n0*m-1-d);
end
N = h^2 * S2 / (delta^2);
N = ceil(N);
N = max(N, n0*ones(1,k));

%%%%%% Step 2
% run Ni-n0 batches for k systems
for i = 1:k
    if N(i) > n0
        Y(:,i,n0+1:N(i)) = kSystem(Xd,i,syspra,N(i)-n0);
    end    
end
Yhat = zeros(m,k);
for i = 1:k
    Yhat(:,i) = mean(Y(:,i,1:N(i)),3);
end

thetahat = zeros(d+1,k);
for i=1:k
    thetahat(:,i) = (Z'*Z)\Z'*Yhat(:,i);  
end

totalN = sum(N) * m;

end