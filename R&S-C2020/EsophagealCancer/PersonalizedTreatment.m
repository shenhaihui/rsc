clear,clc; close all

k = 3; % system number, 1 - no drug, 2 - Aspirin, 3 - Statin
d = 4; % dimension, x = [age,risk,reduction_Aspirin,reduction_Statin]
m = 16; % design points (full factorial design)
Xd = [61 0.1/4 1/4 1/4; ...  % every point, row vector
      61 0.1/4 1/4 3/4; ...
      61 0.1/4 3/4 1/4; ...
      61 0.1/4 3/4 3/4; ...
      61 0.3/4 1/4 1/4; ...
      61 0.3/4 1/4 3/4; ...
      61 0.3/4 3/4 1/4; ...
      61 0.3/4 3/4 3/4; ...
      74 0.1/4 1/4 1/4; ...
      74 0.1/4 1/4 3/4; ...
      74 0.1/4 3/4 1/4; ...
      74 0.1/4 3/4 3/4; ...
      74 0.3/4 1/4 1/4; ...
      74 0.3/4 1/4 3/4; ...
      74 0.3/4 3/4 1/4; ...
      74 0.3/4 3/4 3/4];    
Z = [ones(m,1) Xd];

n0 = 100; % initial sample size
delta = 1/6; % indifferent-zone parameter (2 months)
alpha = 0.05;  % PCS = 1-alpha = 95%

pra.k = k; pra.d = d; pra.m = m; pra.Xd = Xd; pra.Z = Z;
pra.n0 = n0; pra.delta = delta; pra.alpha = alpha;

%%% sove h - EPCS
% h_het = find_h_het_PCSE_EC(pra,1,2)
h_het = 1.9282; % 17145s (285 min), 7 works

%%% carry out procedure FDHet
repn = 300; % the replication number 

% parpool(7); 
% s_pseudo = RandStream('mt19937ar','Seed',100);
% RandStream.setGlobalStream(s_pseudo); 
% 
% %-------- Heter Procedure --------% 
% thetahat_het = zeros(1+d,k,repn);
% N_het = zeros(repn,1);
% tic
% parfor rep = 1:repn
%     [thetahat_het(:,:,rep), N_het(rep)] = select_het_EC(h_het,pra);
% end
% toc
% save('data_het.mat','thetahat_het','N_het');
% fprintf('het is done \n');
% delete(gcp)
load('data_het.mat')


%%%%%  Evaluation  %%%%%%

% "true" value grid
ageR = 55:80;
riskR = 0:0.01:0.1;
AspirinR = 0:0.1:1;
StatinR = 0:0.1:1;
load('QALY_true.mat') % from Brute Force Simulation

% randomly take covariates according to its distribution
s_pseudo = RandStream('mt19937ar','Seed',100);
RandStream.setGlobalStream(s_pseudo);
Xn = 100000;
Xsample = zeros(Xn, 4);

% X1 (Year 2016)
x_prob = [0.058011742, 0.058545686, 0.056257119, 0.055725483, ...
0.054849114, 0.052646386, 0.052030291, 0.049529118, ...
0.047442256, 0.045459377, 0.043798816, 0.042691325, ...
0.041218728, 0.040644411, 0.042182262, 0.030812564, ...
0.030122938, 0.029031999, 0.029415726, 0.025295187, ...
0.022718972, 0.021203217, 0.019726005, 0.018472328, ...
0.016652584, 0.015516365];

x_cumdis = cumsum(x_prob);
for i = 1:Xn
    Xsample(i,1) = find(x_cumdis >= rand, 1) + 54;
end
% X2
Xsample(:,2) = rand(Xn,1)*0.1;
% X3
Xsample(:,3) = TriRnd(0,0.59,1,Xn);
% X4
Xsample(:,4) = TriRnd(0,0.62,1,Xn);

% %%%%%%%%%%%%%%% more specific group %%%%%%%%%%%%%%%
% Xsample(:,3) = 0.9 * ones(Xn,1);
% Xsample(:,4) = 0.2 * ones(Xn,1);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% the "true" best %%%%
% linear interpolation
[G1,G2] = ndgrid(riskR,AspirinR);
GG = [G1(:) G2(:)]; % the same for Statin
y_true = zeros(Xn,3);
for i = 1:Xn  
    X = Xsample(i,:);
    age = X(1);
    
    x1u = ceil(X(2)*100)/100; x1l = floor(X(2)*100)/100;
    x2u = ceil(X(3)*10)/10; x2l = floor(X(3)*10)/10;
    x3u = ceil(X(4)*10)/10; x3l = floor(X(4)*10)/10;    
    if x1u == x1l % risk is integer
        s1 = 0.5;
    else
        s1 = (X(2)-x1u) / (x1l-x1u);
    end
    if x2u == x2l % aspirin is integer
        s2 = 0.5;
    else
        s2 = (X(3)-x2u) / (x2l-x2u);
    end    
    if x3u == x3l % statin is integer
        s3 = 0.5;
    else
        s3 = (X(4)-x3u) / (x3l-x3u);
    end 
    
    % system 1, no drug   
    y_true(i,1) = QALY_nodrug(find(abs(riskR - x1l)<0.001),1,age-54) * s1 + QALY_nodrug(find(abs(riskR - x1u)<0.001),1,age-54) * (1-s1);
    
    % system 2, Aspirin
    row = find( sum(abs(GG-repmat([x1l x2l],size(GG,1),1)),2) < 0.001 );
    y_x1lx2l = QALY_aspirin(row,1,age-54);    
    row = find( sum(abs(GG-repmat([x1u x2l],size(GG,1),1)),2) < 0.001 );
    y_x1ux2l = QALY_aspirin(row,1,age-54);    
    row = find( sum(abs(GG-repmat([x1l x2u],size(GG,1),1)),2) < 0.001 );
    y_x1lx2u = QALY_aspirin(row,1,age-54);       
    row = find( sum(abs(GG-repmat([x1u x2u],size(GG,1),1)),2) < 0.001 );
    y_x1ux2u = QALY_aspirin(row,1,age-54);        
    y_x2l = y_x1lx2l*s1 + y_x1ux2l*(1-s1);
    y_x2u = y_x1lx2u*s1 + y_x1ux2u*(1-s1);    
    y_true(i,2) = y_x2l*s2 + y_x2u*(1-s2);
    
    % system 3, Statin 
    row = find( sum(abs(GG-repmat([x1l x3l],size(GG,1),1)),2) < 0.001 );
    y_x1lx3l = QALY_statin(row,1,age-54);    
    row = find( sum(abs(GG-repmat([x1u x3l],size(GG,1),1)),2) < 0.001 );
    y_x1ux3l = QALY_statin(row,1,age-54);    
    row = find( sum(abs(GG-repmat([x1l x3u],size(GG,1),1)),2) < 0.001 );
    y_x1lx3u = QALY_statin(row,1,age-54);       
    row = find( sum(abs(GG-repmat([x1u x3u],size(GG,1),1)),2) < 0.001 );
    y_x1ux3u = QALY_statin(row,1,age-54);           
    y_x3l = y_x1lx3l*s1 + y_x1ux3l*(1-s1);
    y_x3u = y_x1lx3u*s1 + y_x1ux3u*(1-s1);    
    y_true(i,3) = y_x3l*s3 + y_x3u*(1-s3);        
end
[max_true, choice_true] = max(y_true,[],2);


%-------- Heter Procedure --------% 
CorrR_het_E = zeros(1,repn);
for rep = 1:repn
    %%% evaluate the AEPCS
    yestimate = [ones(Xn,1) Xsample] * thetahat_het(:,:,rep); % for entire test points
    [~, choice] = max(yestimate,[],2);
    index = sub2ind(size(y_true),(1:Xn)',choice);
    correct = choice == choice_true | abs(y_true(index) - max_true) < delta - 1e-6;
    CorrR_het_E(rep) = sum(correct) / Xn;  
end
APCSE_personalized = mean(CorrR_het_E)
% mean(N_het)

%-------- traditional way --------% 
% in expectation sense, which treatment is better? [Always the 3rd one]
% way 1
% mean(y_true)
% 16.9640   17.5274   17.5744

% % way 2
% Eage = 65;
% Erisk = 0.05;
% EAspirin = 0.53;
% EStatin = 0.54;
% X = [Eage, Erisk, EAspirin, EStatin];
% % run simulation for this average patient
% 16.4042   16.9847   17.0233

choice = 3 * ones(Xn,1);
index = sub2ind(size(y_true),(1:Xn)',choice);
correct = choice == choice_true | abs(y_true(index) - max_true) < delta - 1e-6;
APCSE_traditional = sum(correct) / Xn


%%% Regret bar chart under the Selected Treatment Regimen
QALY1 = y_true(:,3);
rep = 1; % arbitrarily choose a replication
yestimate = [ones(Xn,1) Xsample] * thetahat_het(:,:,rep); % for entire test points
[~, choice] = max(yestimate,[],2);
index = sub2ind(size(y_true),(1:Xn)',choice);
QALY2 = y_true(index);

% regret: gap from the true optimal
Regret1 = max_true - QALY1;
Regret2 = max_true - QALY2;

% bar plot/histogram (seperately calculate for 0, relative frequency)
i = 1;
index = Regret1 == 0;
count1(i) = sum(index);
index = Regret2 == 0;
count2(i) = sum(index);
for i = 2:7
    index = (Regret1 > (i-2)*(2/12)) .* (Regret1 <= (i-1)*(2/12));
    count1(i) = sum(index);
    index = (Regret2 > (i-2)*(2/12)) .* (Regret2 <= (i-1)*(2/12));
    count2(i) = sum(index);
end
i = 8;
index = Regret1 > (i-2)*(2/12);
count1(i) = sum(index);
index = Regret2 > (i-2)*(2/12);
count2(i) = sum(index);

count1 = count1 / Xn;
count2 = count2 / Xn;
y = [count1' count2'];
figure, b = bar(y);
b(1).FaceColor = [255,126,13]/255;
b(2).FaceColor = [89,154,212]/255;
set(gca,'xtick',1:8);
xt = {'0','(0,2]','(2,4]','(4,6]','(6,8]','(8,10]','(10,12]','(12,\infty)'};
set(gca,'xticklabel',xt);
xlabel('E[QALYs|X] Regret (unit: month)');
legend('Traditional','Personalized','Location','northeast')