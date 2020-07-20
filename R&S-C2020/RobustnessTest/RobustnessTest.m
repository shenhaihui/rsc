clear, clc

% d = 2,lambda = 0.5, Extreme Design
N_P_R_record = GRF1(0.5, 'extreme');
mean_std_1_1 = [mean(N_P_R_record); std(N_P_R_record)]

% d = 2,lambda = 0.5, Minimax Design
N_P_R_record = GRF1(0.5, 'minimax');
mean_std_1_2 = [mean(N_P_R_record); std(N_P_R_record)]

% d = 2,lambda = 3, Extreme Design
N_P_R_record = GRF1(3, 'extreme');
mean_std_2_1 = [mean(N_P_R_record); std(N_P_R_record)]

% d = 2,lambda = 3, Minimax Design
N_P_R_record = GRF1(3, 'minimax');
mean_std_2_2 = [mean(N_P_R_record); std(N_P_R_record)]

% d = 3,lambda = 0.5, Extreme Design
N_P_R_record = GRF2(0.5, 'extreme');
mean_std_3_1 = [mean(N_P_R_record); std(N_P_R_record)]

% d = 3,lambda = 0.5, Minimax Design
N_P_R_record = GRF2(0.5, 'minimax');
mean_std_3_2 = [mean(N_P_R_record); std(N_P_R_record)]

% d = 3,lambda = 3, Extreme Design
N_P_R_record = GRF2(3, 'extreme');
mean_std_4_1 = [mean(N_P_R_record); std(N_P_R_record)]

% d = 3,lambda = 3, Minimax Design
N_P_R_record = GRF2(3, 'minimax');
mean_std_4_2 = [mean(N_P_R_record); std(N_P_R_record)]