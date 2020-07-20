function N_P_R_record = GRF2(lambda, design)

%% generate surface
x1 = 0:0.01:1;
x2 = 0:0.01:1;
[X1,X2] = ndgrid(x1,x2);
X = [X1(:) X2(:)];
k = size(X,1);

% % parameter of GRF
% tau2 = 1;
% mu = zeros(k,1);
% C = zeros(k,k,2);
% for i = 1:2
%     [Y1,Y2]=meshgrid(X(:,i),X(:,i));
%     C(:,:,i) = Y1- Y2;
% end
% R = tau2 * exp(-lambda * sum(C.^2,3));
% 
% s_pseudo = RandStream('mt19937ar','Seed',0);
% RandStream.setGlobalStream(s_pseudo);
% if lambda == 0.5
%     F_lab = zeros(200,k);
%     for i = 1:200
%         F_lab(i,:) = mvnrnd(mu,R);
%     end
% elseif lambda == 3
%     F_lab = zeros(1500,k);
%     for i = 1:1500
%         F_lab(i,:) = mvnrnd(mu,R);
%     end
% end
% 
% ps = 50; % pool size
% j = 0;
% F = zeros(ps,k);
% for i = 1:1500 %200
%     % f = mvnrnd(mu,R); 
%     f = F_lab(i,:); % use the saved date   
%     % calculate R-squared
%     beta = [ones(k,1) X] \ f';
%     SSE = sum(([ones(k,1) X] * beta - f').^2);
%     SST = sum((mean(f) - f).^2);
%     R2 = 1-SSE/SST;    
%     if R2 >= 0.8          % surface is approximately linear     
%         j = j + 1;
%         F(j,:) = f;
%         if j == ps
%             break;
%         end
%     end
% end

% to save run time, directly use the saved data
if lambda == 0.5
    load('F_05.mat');
elseif lambda == 3
    load('F_3.mat');
end
ps = 50; % pool size

% show the R2
R2_record = zeros(1,50);
for i = 1:50
    f = F(i,:);
    beta = [ones(k,1) X] \ f';
    SSE = sum(([ones(k,1) X] * beta - f').^2);
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
d = 2; % dimension = 2, X1 ~ Unif[0,1], X2 ~ Unif[0,1]
m = 8; % design points
k = 5; % system number

% design matrix
if design == "extreme"
    Xd =[0 0; 0 0; 0 1; 0 1; 1 0; 1 0; 1 1; 1 1];
elseif design == "minimax"
    Xd = [0.1557 0.2086; 0.1557 0.7914; 0.8443 0.2086; 0.8443 0.7914; ...
          0.2468 0.5000; 0.7532 0.5000; 0.5000 0.1794; 0.5000 0.8206];
end

Z = [ones(m,1) Xd];
pra.k = k; pra.d = d; pra.m = m; pra.Xd = Xd; pra.Z = Z;
pra.n0 = n0; pra.delta = delta; pra.alpha = alpha;

%%% PCSE target
% h_hom = find_h_hom_PCSE(pra,1,3)

if design == "extreme"
    h_hom = 1.4718;
elseif design == "minimax"
    h_hom = 2.1468;
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
        index = find( sum(abs(X-repmat(Xsample(j,:),size(X,1),1)),2) < 0.0001 );
        y_true(j,:) = yall(:,index)';
    end
    [max_true, choice_true] = max(y_true,[],2); 

    % the true value at design point
    f_true = zeros(m,k);
    for j = 1:m
        index = find( sum(abs(X-repmat(Xd(j,:),size(X,1),1)),2) < 0.0001 );
        if isempty(index)
            XL1 = Xd(j,1) - mod(Xd(j,1),0.01);
            XR1 = XL1 + 0.01;
            XL2 = Xd(j,2) - mod(Xd(j,2),0.01);
            XR2 = XL2 + 0.01;
            
            indexLL = find( sum(abs(X-repmat([XL1,XL2],size(X,1),1)),2) < 0.0001 );
            indexLR = find( sum(abs(X-repmat([XL1,XR2],size(X,1),1)),2) < 0.0001 );
            fL = 100 * ((XR2 - Xd(j,2)) * yall(:,indexLL)' + (Xd(j,2) - XL2) * yall(:,indexLR)');
            
            indexRL = find( sum(abs(X-repmat([XR1,XL2],size(X,1),1)),2) < 0.0001 );
            indexRR = find( sum(abs(X-repmat([XR1,XR2],size(X,1),1)),2) < 0.0001 );
            fR = 100 * ((XR2 - Xd(j,2)) * yall(:,indexRL)' + (Xd(j,2) - XL2) * yall(:,indexRR)');   
            
            f_true(j,:) = 100 * ((XR1 - Xd(j,1)) * fL + (Xd(j,1) - XL1) * fR);
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