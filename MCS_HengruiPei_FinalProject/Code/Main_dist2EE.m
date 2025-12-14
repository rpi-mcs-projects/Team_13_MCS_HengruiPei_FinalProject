clear all;
User_Position = 0:12.5:300;
Iteration = 8;
Dink_max = 15;
M = 6; K = 4; N = 256;
PBS = db2pow(9); PRIS = db2pow(9);
WBS = db2pow(6); WUE = db2pow(-20); WRIS = db2pow(-20); WRA = db2pow(-20);
xi = 1.1; zeta = 1.1;
sigma1 = 1e-11; % user noise
sigma2 = 1e-11; % RIS noise
fc = 5; % unit: GHz
large_fading_AI = 2; large_fading_DI = 2;
T1 = 1; T2 = 4; T3 = 16; T4 = 64; T5 = 256; % elements per group, L = N/T
C1 = kron(eye(N/T1), ones(T1, 1)); C2 = kron(eye(N/T2), ones(T2, 1));
C3 = kron(eye(N/T3), ones(T3, 1)); C4 = kron(eye(N/T4), ones(T4, 1)); C5 = kron(eye(N/T5), ones(T5, 1)); 
C1_mat = inv(C1'*C1)*C1'; C2_mat = inv(C2'*C2)*C2';
C3_mat = inv(C3'*C3)*C3'; C4_mat = inv(C4'*C4)*C4'; C5_mat = inv(C5'*C5)*C5'; 

R1 = zeros(length(User_Position), Iteration); R2 = R1; R3 = R1; R4 = R1; R5 = R1; R0 = R1;
P1 = zeros(length(User_Position), Iteration); P2 = P1; P3 = P1; P4 = P1; P5 = P1; P0 = P1;

tic;
for u = 1:length(User_Position)
    parfor i = 1:Iteration
        fprintf('index: %d, iter: %d\n', u, i);
        % Position
        [Dis_BStoRIS, Dis_BStoUser, Dis_RIStoUser] = Generate_Position(K, User_Position(u));
        % Channel
        [h_k,f_k,G] = Generate_Channel(K,N,M,large_fading_AI,large_fading_DI,Dis_BStoRIS,Dis_BStoUser,Dis_RIStoUser,fc);
        % Initialization
        W = exp(1j*2*pi*rand(K*M,1))*sqrt(PBS/K/M);
        Theta = diag(exp(1j*2*pi*rand(N,1)));
        % Passive Precoding
        [R0(u, i), P0(u, i)] = Passive_Dinkelbach(xi,M,K,N,PBS+PRIS-N*WRIS,WBS,WUE,WRIS,sigma1,Theta,W,h_k,f_k,G);
        % Active Precoding
        Theta = 1000*Theta;
        [R1(u, i), P1(u, i)] = Active_Dinkelbach(T1,C1,C1_mat,xi,zeta,M,K,N,PBS,PRIS,WBS,WUE,WRIS,WRA,sigma1,sigma2,Theta,W,h_k,f_k,G);
        [R3(u, i), P3(u, i)] = Active_Dinkelbach(T3,C3,C3_mat,xi,zeta,M,K,N,PBS,PRIS,WBS,WUE,WRIS,WRA,sigma1,sigma2,Theta,W,h_k,f_k,G);
    end
end
toc;

EE1 = R1./P1; EE2 = R2./P2; EE3 = R3./P3; EE4 = R4./P4; EE5 = R5./P5; EE0 = R0./P0;
R1 = mean(R1, 2); R2 = mean(R2, 2); R3 = mean(R3, 2); R4 = mean(R4, 2); R5 = mean(R5, 2); R0 = mean(R0, 2);
P1 = mean(P1, 2); P2 = mean(P2, 2); P3 = mean(P3, 2); P4 = mean(P4, 2); P5 = mean(P5, 2); P0 = mean(P0, 2);
EE1 = mean(EE1, 2); EE2 = mean(EE2, 2); EE3 = mean(EE3, 2); EE4 = mean(EE4, 2); EE5 = mean(EE5, 2); EE0 = mean(EE0, 2);

figure;
set(gcf, 'Position', [680,402,560,576]);
subplot(212);
hold on; box on; grid on;
plot(nan, nan, 'LineWidth', 2, 'Marker', '^', 'MarkerSize', 8);
plot(nan, nan, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 8);
plot(nan, nan, 'k--', 'LineWidth', 2);
plot(User_Position, EE3, 'LineWidth', 2, 'Color', '#0072BD');
plot(User_Position, EE1, 'LineWidth', 2, 'Color', '#D95319');
plot(User_Position, EE0, 'k--', 'LineWidth', 2);
scatter(User_Position(1:2:end), EE3(1:2:end), 55, 'MarkerEdgeColor', '#0072BD', 'Marker', '^', 'LineWidth', 2);
scatter(User_Position(1:2:end), EE1(1:2:end), 70, 'MarkerEdgeColor', '#D95319', 'Marker', 'o', 'LineWidth', 2);
a = get(gca, 'xticklabel');
set(gca, 'xticklabel', a, 'FontSize', 11, 'Fontname', 'Times');
xlabel('User position (m)', 'FontSize', 13, 'Fontname', 'Times', 'Interpreter', 'latex');
ylabel('EE (bps/J)', 'FontSize', 13, 'Fontname', 'Times', 'Interpreter', 'latex');
legend('Proposed sub-connected architecture ($T=16$)', 'Fully-connected architecture ($T=1$) [5]', 'Passive RIS [9]', ...
    'FontSize', 12, 'Location', 'southwest', 'Fontname', 'Times', 'Interpreter', 'latex');

subplot(211);
hold on; box on; grid on;
plot(nan, nan, 'LineWidth', 2, 'Marker', '^', 'MarkerSize', 8);
plot(nan, nan, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 8);
plot(nan, nan, 'k--', 'LineWidth', 2);
plot(User_Position, R3, 'LineWidth', 2, 'Color', '#0072BD');
plot(User_Position, R1, 'LineWidth', 2, 'Color', '#D95319');
plot(User_Position, R0, 'k--', 'LineWidth', 2);
scatter(User_Position(1:2:end), R3(1:2:end), 55, 'MarkerEdgeColor', '#0072BD', 'Marker', '^', 'LineWidth', 2);
scatter(User_Position(1:2:end), R1(1:2:end), 70, 'MarkerEdgeColor', '#D95319', 'Marker', 'o', 'LineWidth', 2);
a = get(gca, 'xticklabel');
set(gca, 'xticklabel', a, 'FontSize', 11, 'Fontname', 'Times');
xlabel('User position (m)', 'FontSize', 13, 'Fontname', 'Times', 'Interpreter', 'latex');
ylabel('SE (bps/Hz)', 'FontSize', 13, 'Fontname', 'Times', 'Interpreter', 'latex');
legend('Proposed sub-connected architecture ($T=16$)', 'Fully-connected architecture ($T=1$) [5]', 'Passive RIS [9]', ...
    'FontSize', 12, 'Location', 'southwest', 'Fontname', 'Times', 'Interpreter', 'latex');