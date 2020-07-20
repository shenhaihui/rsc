clear,clc,close all

% fixed parameters
n0 = 50; % initial sample size
delta = 1; % indifferent-zone parameter
alpha = 0.05; % PCS = 1-alpha  

AllCase = cell(2,13);

%%%%%%% Benchmark Case
d = 3; % dimension = 3, X1 ~ Unif[0,1], X2 ~ Unif[0,1], X3 ~ Unif[0,1]
m = 2^d; % design points (full factorial design)
Xd = [ 0  0  0; ...   % every point, row vector
      .5  0  0; ...
       0 .5  0; ...
       0  0 .5; ...
      .5 .5  0; ...
      .5  0 .5; ...
       0 .5 .5; ...
      .5 .5 .5];  
Z = [ones(m,1) Xd];
k = 5; % system number
beta0 = [delta zeros(1,k-1)]; beta = ones(d,k); % system configuration: Slippage Configuration (SC) 
pra.k = k; pra.d = d; pra.m = m; pra.Xd = Xd; pra.Z = Z;
pra.n0 = n0; pra.delta = delta; pra.alpha = alpha;
syspra.beta0 = beta0; syspra.beta = beta;
homo = 1; sigma = ones(1,k) * 10;   % Hom + equal variance (EV)
syspra.homo = homo; syspra.sigma2 = sigma.^2;
%%% sove Xmin
[X1m,X2m,X3m] = ndgrid([0,1],[0,1],[0,1]); % the min is on the corner
Xcorner = [X1m(:) X2m(:) X3m(:)];  % each element, row vector
fcorner = [ones(m,1) Xcorner] / (Z'*Z) * [ones(m,1) Xcorner]';
fcorner = diag(fcorner);
[~, max_l] = max(fcorner);
Xmin = Xcorner(max_l,:);
pra.Xmin = Xmin;
%%% sove h_min
h_hom = find_h_hom_PCSmin(pra,0,20)
h_het = find_h_het_PCSmin(pra,0,20)
pra.h_hom = h_hom;
pra.h_het = h_het;
AllCase{1,1} = pra;
AllCase{2,1} = syspra;


%%%%%%% Case (1) => k = 2
d = 3;
m = 2^d;
Xd = [ 0  0  0; ...
      .5  0  0; ...
       0 .5  0; ...
       0  0 .5; ...
      .5 .5  0; ...
      .5  0 .5; ...
       0 .5 .5; ...
      .5 .5 .5];  
Z = [ones(m,1) Xd];
k = 2;  % <=====
beta0 = [delta zeros(1,k-1)]; beta = ones(d,k); % SC
pra.k = k; pra.d = d; pra.m = m; pra.Xd = Xd; pra.Z = Z;
pra.n0 = n0; pra.delta = delta; pra.alpha = alpha;
syspra.beta0 = beta0; syspra.beta = beta;
homo = 1; sigma = ones(1,k) * 10;   % Hom + EV
syspra.homo = homo; syspra.sigma2 = sigma.^2;
%%% sove Xmin
[X1m,X2m,X3m] = ndgrid([0,1],[0,1],[0,1]); % the min is on the corner
Xcorner = [X1m(:) X2m(:) X3m(:)];  % each element, row vector
fcorner = [ones(m,1) Xcorner] / (Z'*Z) * [ones(m,1) Xcorner]';
fcorner = diag(fcorner);
[~, max_l] = max(fcorner);
Xmin = Xcorner(max_l,:);
pra.Xmin = Xmin;
%%% sove h_min
h_hom = find_h_hom_PCSmin(pra,0,20)
h_het = find_h_het_PCSmin(pra,0,20)
pra.h_hom = h_hom;
pra.h_het = h_het;
AllCase{1,2} = pra;
AllCase{2,2} = syspra;


%%%%%%% Case (2) => k = 8
d = 3;
m = 2^d;
Xd = [ 0  0  0; ...
      .5  0  0; ...
       0 .5  0; ...
       0  0 .5; ...
      .5 .5  0; ...
      .5  0 .5; ...
       0 .5 .5; ...
      .5 .5 .5];  
Z = [ones(m,1) Xd];
k = 8;  % <=====
beta0 = [delta zeros(1,k-1)]; beta = ones(d,k); % SC
pra.k = k; pra.d = d; pra.m = m; pra.Xd = Xd; pra.Z = Z;
pra.n0 = n0; pra.delta = delta; pra.alpha = alpha;
syspra.beta0 = beta0; syspra.beta = beta;
homo = 1; sigma = ones(1,k) * 10;   % Hom + EV
syspra.homo = homo; syspra.sigma2 = sigma.^2;
%%% sove Xmin
[X1m,X2m,X3m] = ndgrid([0,1],[0,1],[0,1]); % the min is on the corner
Xcorner = [X1m(:) X2m(:) X3m(:)];  % each element, row vector
fcorner = [ones(m,1) Xcorner] / (Z'*Z) * [ones(m,1) Xcorner]';
fcorner = diag(fcorner);
[~, max_l] = max(fcorner);
Xmin = Xcorner(max_l,:);
pra.Xmin = Xmin;
%%% sove h_min
h_hom = find_h_hom_PCSmin(pra,0,20)
h_het = find_h_het_PCSmin(pra,0,20)
pra.h_hom = h_hom;
pra.h_het = h_het;
AllCase{1,3} = pra;
AllCase{2,3} = syspra;


%%%%%%% Case (3) => Random generated System
pra = AllCase{1,1};
syspra = AllCase{2,1};
k = 5;
s_pseudo = RandStream('mt19937ar','Seed',12);  % <=====
RandStream.setGlobalStream(s_pseudo); 
beta0 = rand(1,k);
beta = rand(d,k);
syspra.beta0 = beta0; syspra.beta = beta;
AllCase{1,4} = pra;
AllCase{2,4} = syspra;


%%%%%%% Case (4) => Hom + IV
pra = AllCase{1,1};
syspra = AllCase{2,1};
k = 5; 
sigma = 5+((1:k)-1)*10/(k-1);  % Hom + increasing variance (IV) <====
syspra.sigma2 = sigma.^2;
AllCase{1,5} = pra;
AllCase{2,5} = syspra;


%%%%%%% Case (5) => Hom + DV
pra = AllCase{1,1};
syspra = AllCase{2,1};
k = 5;
sigma = 15-((1:k)-1)*10/(k-1); % Hom + decreasing variance (DV) <====
syspra.sigma2 = sigma.^2;
AllCase{1,6} = pra;
AllCase{2,6} = syspra;


%%%%%%% Case (6) => Heteroscedastic Noise
pra = AllCase{1,1};
syspra = AllCase{2,1};
homo = 0;
syspra.homo = homo;
AllCase{1,7} = pra;
AllCase{2,7} = syspra;


%%%%%%% Case (7) => d = 1
d = 1; % dimension = 1, X ~ Unif[0,1]  <====
m = 2^d;
Xd = [0; 1/2];  
Z = [ones(m,1) Xd];
k = 5; % system number
beta0 = [delta zeros(1,k-1)]; beta = ones(d,k); % SC
pra.k = k; pra.d = d; pra.m = m; pra.Xd = Xd; pra.Z = Z;
pra.n0 = n0; pra.delta = delta; pra.alpha = alpha;
syspra.beta0 = beta0; syspra.beta = beta;
homo = 1; sigma = ones(1,k) * 10;   % Hom + EV
syspra.homo = homo; syspra.sigma2 = sigma.^2;
%%% sove Xmin
[X1m] = ndgrid([0,1]); % the min is on the corner
Xcorner = X1m(:);  % each element, row vector
fcorner = [ones(m,1) Xcorner] / (Z'*Z) * [ones(m,1) Xcorner]';
fcorner = diag(fcorner);
[~, max_l] = max(fcorner);
Xmin = Xcorner(max_l,:);
pra.Xmin = Xmin;
%%% sove h_min
h_hom = find_h_hom_PCSmin(pra,0,20)
h_het = find_h_het_PCSmin(pra,0,20)
pra.h_hom = h_hom;
pra.h_het = h_het;
AllCase{1,8} = pra;
AllCase{2,8} = syspra;


%%%%%%% Case (8) => d = 5
d = 5; % dimension = 3, X1 ~ Unif[0,1], X2 ~ Unif[0,1], X3 ~ Unif[0,1]  <====
% X4 ~ Unif[0,1], X5 ~ Unif[0,1]
m = 2^d; % design points (full factorial design)
[X1m,X2m,X3m,X4m,X5m] = ndgrid([0,0.5],[0,0.5],[0,0.5],[0,0.5],[0,0.5]);
Xd = [X1m(:) X2m(:) X3m(:) X4m(:) X5m(:)];
Z = [ones(m,1) Xd];
k = 5; % system number
beta0 = [delta zeros(1,k-1)]; beta = ones(d,k); % SC
pra.k = k; pra.d = d; pra.m = m; pra.Xd = Xd; pra.Z = Z;
pra.n0 = n0; pra.delta = delta; pra.alpha = alpha;
syspra.beta0 = beta0; syspra.beta = beta;
homo = 1; sigma = ones(1,k) * 10;   % Hom + EV
syspra.homo = homo; syspra.sigma2 = sigma.^2;
%%% sove Xmin
[X1m,X2m,X3m,X4m,X5m] = ndgrid([0,1],[0,1],[0,1],[0,1],[0,1]); % the min is on the corner
Xcorner = [X1m(:) X2m(:) X3m(:) X4m(:) X5m(:)];  % each element, row vector
fcorner = [ones(m,1) Xcorner] / (Z'*Z) * [ones(m,1) Xcorner]';
fcorner = diag(fcorner);
[~, max_l] = max(fcorner);
Xmin = Xcorner(max_l,:);
pra.Xmin = Xmin;
%%% sove h_min
h_hom = find_h_hom_PCSmin(pra,0,20)
h_het = find_h_het_PCSmin(pra,0,20)
pra.h_hom = h_hom;
pra.h_het = h_het;
AllCase{1,9} = pra;
AllCase{2,9} = syspra;


%%%%%%% Case (9) => X_i follows truncated normal distribution with corr=0.5
d = 3;
m = 2^d; % design points (full factorial design)
Xd = [ 0  0  0; ...   % every point, row vector
      .5  0  0; ...
       0 .5  0; ...
       0  0 .5; ...
      .5 .5  0; ...
      .5  0 .5; ...
       0 .5 .5; ...
      .5 .5 .5];  
Z = [ones(m,1) Xd];
k = 5; % system number
beta0 = [delta zeros(1,k-1)]; beta = ones(d,k); %SC
pra.k = k; pra.d = d; pra.m = m; pra.Xd = Xd; pra.Z = Z;
pra.n0 = n0; pra.delta = delta; pra.alpha = alpha;
syspra.beta0 = beta0; syspra.beta = beta;
homo = 1; sigma = ones(1,k) * 10;   % Hom + equal variance (EV)
syspra.homo = homo; syspra.sigma2 = sigma.^2;
%%% sove Xmin
[X1m,X2m,X3m] = ndgrid([0,1],[0,1],[0,1]); % the min is on the corner
Xcorner = [X1m(:) X2m(:) X3m(:)];  % each element, row vector
fcorner = [ones(m,1) Xcorner] / (Z'*Z) * [ones(m,1) Xcorner]';
fcorner = diag(fcorner);
[~, max_l] = max(fcorner);
Xmin = Xcorner(max_l,:);
pra.Xmin = Xmin;
%%% sove h_min
h_hom = find_h_hom_PCSmin(pra,0,20)
h_het = find_h_het_PCSmin(pra,0,20)
pra.h_hom = h_hom;
pra.h_het = h_het;
AllCase{1,10} = pra;
AllCase{2,10} = syspra;


%%%%%%% large-scale problems
%%%%%%% Case (10) => k = 100
d = 3; % dimension = 3, X1 ~ Unif[0,1], X2 ~ Unif[0,1], X3 ~ Unif[0,1]
m = 2*d; % design points (Latin hypercube sample)
s_pseudo = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s_pseudo);    
Xd = lhsdesign(m,d);
Z = [ones(m,1) Xd];
k = 100; % system number       <====
beta0 = [delta zeros(1,k-1)]; beta = ones(d,k); % SC
pra.k = k; pra.d = d; pra.m = m; pra.Xd = Xd; pra.Z = Z;
pra.n0 = n0; pra.delta = delta; pra.alpha = alpha;
syspra.beta0 = beta0; syspra.beta = beta;
homo = 1; sigma = ones(1,k) * 10;   % Hom + EV
syspra.homo = homo; syspra.sigma2 = sigma.^2;
%%% sove Xmin
lb = [0 0 0];
ub = [1 1 1];
fun = @(x) -[1 x] / (Z'*Z) * [1 x]';
Xmin = fmincon(fun,0.5*ones(1,d),[],[],[],[],lb,ub);
pra.Xmin = Xmin;
%%% sove h_min
h_hom = find_h_hom_PCSmin(pra,0,20)
h_het = find_h_het_PCSmin(pra,0,20)
pra.h_hom = h_hom;
pra.h_het = h_het;
AllCase{1,11} = pra;
AllCase{2,11} = syspra;


%%%%%%% Case (11) => d = 49+1
d = 49; % dimension = 49, X_i ~ Unif[0,1]  <====
m = 2*d; % design points (Latin hypercube sample)
s_pseudo = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s_pseudo);    
Xd = lhsdesign(m,d);
Z = [ones(m,1) Xd];
k = 5; % system number
beta0 = [delta zeros(1,k-1)]; beta = ones(d,k); % SC
pra.k = k; pra.d = d; pra.m = m; pra.Xd = Xd; pra.Z = Z;
pra.n0 = n0; pra.delta = delta; pra.alpha = alpha;
syspra.beta0 = beta0; syspra.beta = beta;
homo = 1; sigma = ones(1,k) * 10;   % Hom + EV
syspra.homo = homo; syspra.sigma2 = sigma.^2;
%%% sove Xmin
lb = zeros(1,d);
ub = ones(1,d);
fun = @(x) -[1 x] / (Z'*Z) * [1 x]';
Xmin = fmincon(fun,0.5*ones(1,d),[],[],[],[],lb,ub);
pra.Xmin = Xmin;
%%% sove h_min
h_hom = find_h_hom_PCSmin(pra,0,20)
h_het = find_h_het_PCSmin(pra,0,20)
pra.h_hom = h_hom;
pra.h_het = h_het;
AllCase{1,12} = pra;
AllCase{2,12} = syspra;


%%%%%%% Case (12) => d = 49+1, k = 100
d = 49; % dimension = 49, X_i ~ Unif[0,1]  <====
m = 2*d; % design points (Latin hypercube sample)
s_pseudo = RandStream('mt19937ar','Seed',0);
RandStream.setGlobalStream(s_pseudo);    
Xd = lhsdesign(m,d);
Z = [ones(m,1) Xd];
k = 100; % system number
beta0 = [delta zeros(1,k-1)]; beta = ones(d,k); % SC
pra.k = k; pra.d = d; pra.m = m; pra.Xd = Xd; pra.Z = Z;
pra.n0 = n0; pra.delta = delta; pra.alpha = alpha;
syspra.beta0 = beta0; syspra.beta = beta;
homo = 1; sigma = ones(1,k) * 10;   % Hom + EV
syspra.homo = homo; syspra.sigma2 = sigma.^2;
%%% sove Xmin
lb = zeros(1,d);
ub = ones(1,d);
fun = @(x) -[1 x] / (Z'*Z) * [1 x]';
Xmin = fmincon(fun,0.5*ones(1,d),[],[],[],[],lb,ub);
pra.Xmin = Xmin;
%%% sove h_min
h_hom = find_h_hom_PCSmin(pra,0,20)
h_het = find_h_het_PCSmin(pra,0,20)
pra.h_hom = h_hom;
pra.h_het = h_het;
AllCase{1,13} = pra;
AllCase{2,13} = syspra;


%%%%%%% Start to run all problems
CorrR_hom_all = zeros(13,3);
N_hom_all = zeros(13,1);
CorrR_het_all = zeros(13,3);
N_het_all = zeros(13,1);

for i = 1:13
    % control random seed
    s_pseudo = RandStream('mt19937ar','Seed',100);
    RandStream.setGlobalStream(s_pseudo); 
    
    pra = AllCase{1,i};
    syspra = AllCase{2,i};
    h_hom = pra.h_hom;
    h_het = pra.h_het;
    d = pra.d;
    k = pra.k;
    Xmin = pra.Xmin;
    
    % verifying samples
    Xn = 100000;
    if i ~= 9+1 % X_i ~ Unif[0,1]
        Xsample = rand(Xn,pra.d);  
    else % X_i ~ truncated normal distribution with corr=0.5
        mu = [0.5;0.5;0.5];
        Sigma = [1 0.5 0.5;0.5 1 0.5; 0.5 0.5 1];
        Xsample = [];
        while size(Xsample,1) < Xn
            Xsample = [Xsample; mvnrnd(mu,Sigma,Xn)];
            index = ones(size(Xsample,1),1);
            for j = 1:3
                index = index .* (Xsample(:,j) <= 1) .* (Xsample(:,j) >= 0);
            end
            Xsample = Xsample(find(index == 1),:);        
        end
        Xsample = Xsample(1:Xn,:);        
    end
    
    if i == 3+1  % general configuration
        beta0 = syspra.beta0;
        beta = syspra.beta;
        y_true = [ones(Xn,1) Xsample] * [beta0; beta];
        [max_true, choice_true] = max(y_true,[],2); 
        
        y_true1 = [1 Xmin] * [beta0; beta];
        [max_true1, choice_true1] = max(y_true1,[],2);
    end
    
    repn = 10000; % the replication number    
    
    %-------- Homo Procedure --------%
    thetahat_hom = zeros(1+d,k,repn);
    N_hom = zeros(repn,1);
    CorrR_hom_E = zeros(1,repn);
    CorrR_hom_min = zeros(1,repn);
    CorrR_hom_min_new = zeros(Xn,1);
    tStart = tic;
    tSim = 0;
    for rep = 1:repn
        tic
        [thetahat_hom(:,:,rep), N_hom(rep)] = select_hom(h_hom,pra,syspra);
        tSim = tSim + toc;
             
        %%% evaluate the PCSE
        yestimate = [ones(Xn,1) Xsample] * thetahat_hom(:,:,rep); % for entire test points
        [~, choice] = max(yestimate,[],2);  
        if i ~= 3+1  % SC case, system 1 is the unique best
            correct = choice == 1;
        else         % general configuration
            index = sub2ind(size(y_true),(1:Xn)',choice);
            correct = choice == choice_true | abs(y_true(index) - max_true) < delta - 1e-6;                
        end
        CorrR_hom_E(rep) = sum(correct) / Xn;             
        
        %%% evaluate the PCSmin
        yestimate = [1 Xmin] * thetahat_hom(:,:,rep); % for the worst point
        [~, choice] = max(yestimate,[],2);        
        if i ~= 3+1  % SC case, system 1 is the unique best
            CorrR_hom_min(rep) = choice == 1;
        else         % general configuration
            CorrR_hom_min(rep) = choice == choice_true1 | abs(y_true1(choice) - max_true1) < delta - 1e-6;  
        end        
        
        %%% evaluate the new PCSmin
        CorrR_hom_min_new = CorrR_hom_min_new + correct;        
    end
    tEnd = toc(tStart);
    fprintf('Case %d, Procedure TS, simulation run time %d seconds, total run time %d seconds \n',i-1,tSim,tEnd);
    CorrR_hom_all(i,:) = [mean(CorrR_hom_E) mean(CorrR_hom_min) min(CorrR_hom_min_new/repn)];
    N_hom_all(i) = mean(N_hom);    
    
    
    %-------- Heter Procedure --------% 
    thetahat_het = zeros(1+d,k,repn);
    N_het = zeros(repn,1);
    CorrR_het_E = zeros(1,repn);
    CorrR_het_min = zeros(1,repn);
    CorrR_het_min_new = zeros(Xn,1);    
    tStart = tic;
    tSim = 0;    
    for rep = 1:repn
        tic
        [thetahat_het(:,:,rep), N_het(rep)] = select_het(h_het,pra,syspra);
        tSim = tSim + toc;
                
        %%% evaluate the PCSE
        yestimate = [ones(Xn,1) Xsample] * thetahat_het(:,:,rep); % for entire test points
        [~, choice] = max(yestimate,[],2);
        if i ~= 3+1  % SC case, system 1 is the unique best
            correct = choice == 1;
        else         % general configuration
            index = sub2ind(size(y_true),(1:Xn)',choice);
            correct = choice == choice_true | abs(y_true(index) - max_true) < delta - 1e-6;
        end
        CorrR_het_E(rep) = sum(correct) / Xn;  
        
        %%% evaluate the PCSmin
        yestimate = [1 Xmin] * thetahat_het(:,:,rep); % for the worst points
        [~, choice] = max(yestimate,[],2);        
        if i ~= 3+1  % SC case, system 1 is the unique best
            CorrR_het_min(rep) = choice == 1;
        else         % general configuration
            CorrR_het_min(rep) = choice == choice_true1 | abs(y_true1(choice) - max_true1) < delta - 1e-6;  
        end        
        
        %%% evaluate the new PCSmin
        CorrR_het_min_new = CorrR_het_min_new + correct;          
    end
    tEnd = toc(tStart);
    fprintf('Case %d, Procedure TS+, simulation run time %d seconds, total run time %d seconds \n',i-1,tSim,tEnd);
    CorrR_het_all(i,:) = [mean(CorrR_het_E) mean(CorrR_het_min) min(CorrR_het_min_new/repn)];
    N_het_all(i) = mean(N_het);        
end