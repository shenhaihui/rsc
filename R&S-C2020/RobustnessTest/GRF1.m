function N_P_R_record = GRF1(lambda, design)

%% generate true surfaces
x = 0:0.01:1;
k = length(x);

% parameter of GRF
tau2 = 1;
mu = zeros(k,1);
[X,Y]=meshgrid(x,x);
C = X - Y;
R = tau2 * exp(-lambda*C.^2);

s_pseudo = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s_pseudo); 
ps = 50; % pool size
j = 0;
F = zeros(ps,k);
for i = 1:500
    f = mvnrnd(mu,R);
    Corr = corrcoef(x,f);
    if abs(Corr(1,2)) >= sqrt(0.8)  % surface is approximately linear
        j = j + 1;
        F(j,:) = f;
        if j == ps
            break;
        end
    end
end

% show the R2
R2_record = zeros(1,50);
for i = 1:50
    f = F(i,:);
    beta = [ones(k,1) x'] \ f';
    SSE = sum(([ones(k,1) x'] * beta - f').^2);
    SST = sum((mean(f) - f).^2);
    R2 = 1-SSE/SST;
    R2_record(i) = R2;
end
[mean(R2_record) var(R2_record)]


%% R&S-C
% fixed parameters
n0 = 50; % initial sample size
delta = 0.2; % indifferent-zone parameter
alpha = 0.05; % PCS = 1-alpha  
d = 1; % dimension = 1, X1 ~ Unif[0,1]
m = 4; % design points
k = 5; % system number

% design matrix
if design == "extreme"
    Xd = [0; 0; 1; 1];
elseif design == "minimax"
    Xd = [1/8; 3/8; 5/8; 7/8];
end

Z = [ones(m,1) Xd];
pra.k = k; pra.d = d; pra.m = m; pra.Xd = Xd; pra.Z = Z;
pra.n0 = n0; pra.delta = delta; pra.alpha = alpha;

%%% PCSE target
% h_hom = find_h_hom_PCSE(pra,1,4)

if design == "extreme"
    h_hom = 1.8543;
elseif design == "minimax"
    h_hom = 2.3902;
end

pra.h_hom = h_hom;


%%
N_P_R_record = zeros(100,3);
for i = 1:100
    % control random seed
    s_pseudo = RandStream('mt19937ar','Seed',i);
    RandStream.setGlobalStream(s_pseudo);    
    
    % generate a problem k = 5
    index = zeros(k,1);
    for j = 1:k
        index2 = unidrnd(ps);
        while sum(index == index2) > 0
            index2 = unidrnd(ps);
        end
        index(j) = index2;
    end
    yall = F(index,:);
    
    % verifying samples
    Xn = 10000;
    Xsample = (unidrnd(101,Xn,d)-1) * 0.01;  % range: 0:0.01:1
    
    % calculate the true best
    y_true = zeros(Xn,k);
    for j = 1:Xn
        index = find(abs(x - Xsample(j)) < 0.0001);
        y_true(j,:) = yall(:,index)';
    end
    [max_true, choice_true] = max(y_true,[],2); 

    % the true value at design point
    f_true = zeros(m,k);
    for j = 1:m
        index = find(abs(x - Xd(j)) < 0.0001);
        if isempty(index)
            XL = Xd(j) - mod(Xd(j),0.01);
            XR = XL + 0.01;
            indexL = find(abs(x - XL) < 0.0001);
            indexR = find(abs(x - XR) < 0.0001);
            f_true(j,:) = 100 * ((XR - Xd(j)) * yall(:,indexL)' + (Xd(j) - XL) * yall(:,indexR)');
        else
            f_true(j,:) = yall(:,index)';
        end
    end
    
    repn = 1000; % the replication number    
    
    %-------- Homo Procedure --------%
    thetahat_hom = zeros(1+d,k,repn);
    N_hom = zeros(repn,1);
    CorrR_hom_E = zeros(1,repn);
    Regret = zeros(repn,1);

    for rep = 1:repn
        [thetahat_hom(:,:,rep), N_hom(rep)] = select_hom(h_hom,pra,f_true);
             
        %%% evaluate the PCSE
        yestimate = [ones(Xn,1) Xsample] * thetahat_hom(:,:,rep); % for entire test points
        [~, choice] = max(yestimate,[],2);  

        index = sub2ind(size(y_true),(1:Xn)',choice);
        correct = choice == choice_true | abs(y_true(index) - max_true) < delta - 1e-6;
        CorrR_hom_E(rep) = sum(correct) / Xn;             
        
        %%% evaluate regret
        Regret(rep) = mean(max_true - y_true(index));

    end
    N_P_R_record(i,:) = [mean(N_hom)  mean(CorrR_hom_E)  mean(Regret)];    
end

end