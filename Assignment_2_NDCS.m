clear all
close all
clc

%
ID = [5 3 0 8 7 3 9];
a = ID(1);
b = ID(3);
c = ID(end);
A = [0.3+a-b 0.5-c; 0 1];
B = [0; 1];

sys_c1 = ss(A,B,[],[]);
Poles = [-1+2*i -1-2*i];
K_1 = place(sys_c1.A, sys_c1.B, Poles);

A2 = (1/3) * A;
B2 = sys_c1.B;
sys_c2 = ss(A2,B2,[],[]);
K_2 = K_1 * 0.4;
h_stable_0 = [];
h_stable_h = [];

h = 0.01;

for h = 0.0005:0.0005:0.4
    h
    F_1 = expm(sys_c1.A * h);
    G_1 = inv(sys_c1.A) * ( expm(sys_c1.A * h) - eye(2)) * sys_c1.B;
    
    F_2 = expm(A2 * h);
    G_2 = inv(A2) * ( expm(A2 * h) - eye(2)) * B2;
    
    % To hold 
    F_cl_1_h_sent = [F_1-G_1*K_1 zeros(2,1); 
                    -K_1          0];
    
    F_cl_1_h_lost = [F_1     G_1; 
                     0    0    1];
    
    F_cl_2_h_sent = [F_2-G_2*K_2 zeros(2,1); 
                    -K_2            0];
    
    F_cl_2_h_lost = [F_2     G_2; 
                     0    0    1];
    % To zero
    F_cl_1_0_sent = F_1-G_1*K_1; 
    F_cl_2_0_sent = F_2-G_2*K_2; 
    F_cl_1_0_lost = blkdiag(F_1);
    F_cl_2_0_lost = blkdiag(F_2);
    
    % Now there are two possible sequences for system 1 [Send, Send, Send] or [Send, Send, Lose];
    % Now there is one sequence for system 2: [send lose]
    Sys1_seq1_h = F_cl_1_h_sent^3;
    Sys1_seq1_0 = F_cl_1_0_sent^3;
    Sys1_seq2_h = F_cl_1_h_lost * F_cl_1_h_sent^2;
    Sys1_seq2_0 = F_cl_1_0_lost * F_cl_1_0_sent^2;

    Sys2_0 = F_cl_2_0_lost * F_cl_2_0_sent;
    Sys2_h = F_cl_2_h_lost * F_cl_2_h_sent;
 
    % Now perform stability analysis by looking at LMI = A'PA - P <= -Q;
    % TO ZERO
    P_0 = sdpvar(2,2);
    obj_0 = 0;
    Q0 = eye(2);
    cons_0 = [P_0 >= eye(2)*0.01;
            Sys1_seq1_0' * P_0 * Sys1_seq1_0 - P_0 <= - Q0;
            Sys1_seq2_0' * P_0 * Sys1_seq2_0 - P_0 <= - Q0;
            Sys2_0' * P_0 * Sys2_0 - P_0 <= - Q0;
        ];

     S = sdpsettings('solver', 'sedumi');
     result_0 = optimize(cons_0,obj_0, S);
     disp(result_0.info)
     value(P_0);

     if result_0.problem == 0
        h_stable_0(end+1) = 1;
     else
        h_stable_0(end+1) = 0;
     end
    
     % TO HOLD
    P_h = sdpvar(3,3);
    obj_h = 0;
    Qh = eye(3);
    cons_h = [P_h >= eye(3)*0.01;
            Sys1_seq1_h' * P_h * Sys1_seq1_h - P_h <= - Qh;
            Sys1_seq2_h' * P_h * Sys1_seq2_h - P_h <= - Qh;
            Sys2_h' * P_h * Sys2_h - P_h <= - Qh;
        ]
     S = sdpsettings('solver', 'sedumi');
     result_h = optimize(cons_h,obj_h, S);
     disp(result_h.info)
     value(P_h);
     

     if result_h.problem == 0
        h_stable_h(end+1) = 1;
     else
        h_stable_h(end+1) = 0;
     end
end
%%
time = 0.0005:0.0005:h;
figure(1)
plot(time, h_stable_0, 'r')
hold on 
scatter(time, h_stable_0, 'r')
legend('To-zero', 'Interpreter', 'latex')
grid on 

figure(2)
scatter(time, h_stable_h, 'b')
hold on
plot(time, h_stable_h, 'b')
legend('To-hold', 'Interpreter', 'latex')
grid on

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h_stable_02 = [];
h_stable_h2 = [];
% LMIs
P = [0  0.99    0   0    0  0.01    0;
    0   0     0.99 0   0.01     0   0  ;
    0.49  0 0  0.02 0 0 0.49;
    0 0.99 0 0 0 0.01 0;
    0.49 0 0 0.02 0 0 0.49;
    0  0 0.99 0 0.01 0 0;
    0 0.99 0 0 0 0.01  0; 
    ];

for h = 0.0005:0.0005:0.8
    % System definition
    h
    F_1 = expm(sys_c1.A * h);
    G_1 = inv(sys_c1.A) * ( expm(sys_c1.A * h) - eye(2)) * sys_c1.B;
    F_2 = expm(A2 * h);
    G_2 = inv(A2) * ( expm(A2 * h) - eye(2)) * B2;
    
    % To hold 
    F_cl_1_h_sent = [F_1-G_1*K_1 zeros(2,1); 
                    -K_1          0];
    F_cl_1_h_lost = [F_1     G_1; 
                     0    0    1];
    F_cl_2_h_sent = [F_2-G_2*K_2 zeros(2,1); 
                    -K_2            0];
    F_cl_2_h_lost = [F_2     G_2; 
                     0    0    1];
  
    % Now perform stability analysis by looking at LMI 
    % TO HOLD
    P1 = sdpvar(6,6);
    P2 = sdpvar(6,6);
    P3 = sdpvar(6,6);
    P4 = sdpvar(6,6);
    P5 = sdpvar(6,6);
    P6 = sdpvar(6,6);
    P7 = sdpvar(6,6);

    A_1 = blkdiag(F_cl_1_h_sent, F_cl_2_h_lost);
    A_2 = blkdiag(F_cl_1_h_sent, eye(3)*0);
    A_3 = A_2;
    A_4 = blkdiag(F_cl_1_h_lost, F_cl_2_h_lost);
    A_5 = blkdiag(F_cl_1_h_lost, eye(3)*0);
    A_6 = A_5;
    A_7 = blkdiag(F_cl_1_h_lost, F_cl_2_h_sent);

    % Optimization 
    Qh = eye(6)*0.001;

    cons_hold = [
                P1 >= eye(6)*0.001;
                P2 >= eye(6)*0.001;
                P3 >= eye(6)*0.001;
                P4 >= eye(6)*0.001;
                P5 >= eye(6)*0.001;
                P6 >= eye(6)*0.001;
                P7 >= eye(6)*0.001;
                
                P1 - A_1'*(P(1,1)*P1 + P(1,2) * P2 + P(1,3) * P3 + P(1,4) * P4 + P(1,5) * P5 + P(1,6) * P6 + P(1,7) * P7) * A_1 >= Qh;
                P2 - A_2'*(P(2,1)*P1 + P(2,2) * P2 + P(2,3) * P3 + P(2,4) * P4 + P(2,5) * P5 + P(2,6) * P6 + P(2,7) * P7) * A_2 >= Qh;
                P3 - A_3'*(P(3,1)*P1 + P(3,2) * P2 + P(3,3) * P3 + P(3,4) * P4 + P(3,5) * P5 + P(3,6) * P6 + P(3,7) * P7) * A_3 >= Qh;
                P4 - A_4'*(P(4,1)*P1 + P(4,2) * P2 + P(4,3) * P3 + P(4,4) * P4 + P(4,5) * P5 + P(4,6) * P6 + P(4,7) * P7) * A_4 >= Qh;
                P5 - A_5'*(P(5,1)*P1 + P(5,2) * P2 + P(5,3) * P3 + P(5,4) * P4 + P(5,5) * P5 + P(5,6) * P6 + P(5,7) * P7) * A_5 >= Qh;
                P6 - A_6'*(P(6,1)*P1 + P(6,2) * P2 + P(6,3) * P3 + P(6,4) * P4 + P(6,5) * P5 + P(6,6) * P6 + P(7,7) * P7) * A_6 >= Qh;
                P7 - A_7'*(P(7,1)*P1 + P(7,2) * P2 + P(7,3) * P3 + P(7,4) * P4 + P(7,5) * P5 + P(7,6) * P6 + P(7,7) * P7) * A_7 >= Qh;
                ];

     S = sdpsettings('solver', 'sedumi');
     result_h = optimize(cons_hold, obj_0, S);
     disp(result_h.info)

     if result_h.problem == 0
        h_stable_h2(end+1) = 1;
     else
        h_stable_h2(end+1) = 0;
     end
end
%%

time2 = 0.0005:0.0005:h;
figure(1)
plot(time2, h_stable_h2)
legend('To-hold', 'Interpreter', 'latex')
grid on

%% Question 3

a = ID(1);
b = ID(3);
c = ID(end);
A = [0.3+a-b 0.5-c; 0 1];
B = [0; 1];

[Q, J] = eig(A);
lambda1 = J(1,1);
lambda2 = J(2,2);
% expm(A*h)
% J = jordan(A)
F13 = [];
F23 = [];
tau = [];
alpha1_v = [];
alpha2_v = [];

h = 0.25;
for t=0:0.005:h-0.005
    tau(end+1) = t;
    % F matrix
    F0 = [expm(A*h) zeros(2,1); 0 0 0];
    alpha1 = exp(lambda1*(h-t));
    F_1 = inv(A) * [(-1 + exp(lambda1*t)) 0; 0 0] * B;
    F1 = [eye(2)*0 F_1; 0 0 0];
    alpha2 = exp(-lambda2 * (h - t));
    F_2 = inv(A) * [0 0; 0 (-1+exp(lambda2*t))] * B;
    F2 = [eye(2)*0 F_2; 0 0 0];
    F = F0 + alpha1*F1 + alpha2*F2;

    F13(end+1) = F(1,3);
    F23(end+1) = F(2,3);
    % G matrix
    G0 = [0; 0; 1];
    G_1 = inv(A) * [(1-exp(-lambda1*(h-t))) 0; 0 0] * B;
    G1 = [G_1; 0];
    G_2 = inv(A) * [0 0; 0 (1-exp(-lambda2*(h-t)))] * B;
    G2 = [G_2; 0];
    G = G0 + alpha1*G1 + alpha2 * G2;
    alpha2_v(end+1) = alpha2;
    alpha1_v(end+1) = alpha1;

end
alpha1_max = max(alpha1_v)
alpha1_min = min(alpha1_v)
alpha2_max = max(alpha2_v)
alpha2_min = min(alpha2_v)
X = [alpha1_min alpha1_max, alpha2_max alpha2_min]
Y = [alpha1_min alpha1_max, alpha2_max alpha2_min]
figure(31)
plot(F13,F23)
hold on
area(X,Y)
% plot(F23)
