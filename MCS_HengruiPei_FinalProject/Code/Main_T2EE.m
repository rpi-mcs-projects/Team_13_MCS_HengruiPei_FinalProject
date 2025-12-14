clear all;
User_Position = 100;
Iteration = 16;
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

R1 = zeros(Iteration, 1); R2 = R1; R3 = R1; R4 = R1; R5 = R1; R0 = R1;
P1 = zeros(Iteration, 1); P2 = P1; P3 = P1; P4 = P1; P5 = P1; P0 = P1;

tic;
% Position
[Dis_BStoRIS, Dis_BStoUser, Dis_RIStoUser] = Generate_Position(K, User_Position);
parfor i = 1:Iteration
    fprintf('iter: %d\n', i);
    % Channel
    [h_k,f_k,G] = Generate_Channel(K,N,M,large_fading_AI,large_fading_DI,Dis_BStoRIS,Dis_BStoUser,Dis_RIStoUser,fc);
    % Initialization
    W = exp(1j*2*pi*rand(K*M,1))*sqrt(PBS/K/M);
    Theta = diag(exp(1j*2*pi*rand(N,1)));
    Theta = 1000*Theta;
    % Precoding
    [R1(i), P1(i)] = Active_Dinkelbach(T1,C1,C1_mat,xi,zeta,M,K,N,PBS,PRIS,WBS,WUE,WRIS,WRA,sigma1,sigma2,Theta,W,h_k,f_k,G);
    [R2(i), P2(i)] = Active_Dinkelbach(T2,C2,C2_mat,xi,zeta,M,K,N,PBS,PRIS,WBS,WUE,WRIS,WRA,sigma1,sigma2,Theta,W,h_k,f_k,G);
    [R3(i), P3(i)] = Active_Dinkelbach(T3,C3,C3_mat,xi,zeta,M,K,N,PBS,PRIS,WBS,WUE,WRIS,WRA,sigma1,sigma2,Theta,W,h_k,f_k,G);
    [R4(i), P4(i)] = Active_Dinkelbach(T4,C4,C4_mat,xi,zeta,M,K,N,PBS,PRIS,WBS,WUE,WRIS,WRA,sigma1,sigma2,Theta,W,h_k,f_k,G);
    [R5(i), P5(i)] = Active_Dinkelbach(T5,C5,C5_mat,xi,zeta,M,K,N,PBS,PRIS,WBS,WUE,WRIS,WRA,sigma1,sigma2,Theta,W,h_k,f_k,G);
end
toc;

EE1 = R1./P1; EE2 = R2./P2; EE3 = R3./P3; EE4 = R4./P4; EE5 = R5./P5;
EE1 = mean(EE1); EE2 = mean(EE2); EE3 = mean(EE3); EE4 = mean(EE4); EE5 = mean(EE5);
figure;
hold on; box on; grid on;
plot(1:5, [EE1, EE2, EE3, EE4, EE5], 'LineWidth', 2, 'Marker', 's', 'MarkerSize', 8);
a = get(gca, 'xticklabel');
set(gca, 'xticklabel', a, 'FontSize', 11, 'Fontname', 'Times');
xticks([1 2 3 4 5]);
xticklabels({num2str(T1), num2str(T2), num2str(T3), num2str(T4), num2str(T5)});
xlabel('$T$', 'FontSize', 13, 'Fontname', 'Times', 'Interpreter', 'latex');
ylabel('EE (bit/Joule)', 'FontSize', 13, 'Fontname', 'Times', 'Interpreter', 'latex');
text(1, EE1, num2str(EE1), 'FontSize', 15, 'FontName', 'times');
text(2, EE2, num2str(EE2), 'FontSize', 15, 'FontName', 'times');
text(3, EE3, num2str(EE3), 'FontSize', 15, 'FontName', 'times');
text(4, EE4, num2str(EE4), 'FontSize', 15, 'FontName', 'times');
text(4.8, EE5, num2str(EE5), 'FontSize', 15, 'FontName', 'times');