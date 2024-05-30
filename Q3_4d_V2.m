clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Question 3a
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ID = [5 3 0 8 7 3 9];
% ID = [5 2 7 4 2 3 0]; % kerims student number
a = ID(1);
b = ID(3);
c = ID(end);
A = [0.3+a-b 0.5-c; 0 1];
B = [0; 1];

sys_c1 = ss(A,B,[],[]);
Poles = [-1+2*i -1-2*i];
K_1 = place(sys_c1.A, sys_c1.B, Poles);

h_stable_0 = [];
h_stable_h = [];

[V,D] = eig(A);
% S = V
inv_A = A^(-1);
F23 = [];
F13 = [];
alpha1_v = [];
alpha2_v = [];
F = [];
G = [];
K_cl2 = [K_1 0];
K2 = [];
h_index = 0;
tau_index = 0;

for h = 0.005: 0.005:0.7;
    h_index = h_index +1;
    tau_index = 0;
    h
    H(h_index) = h;
    for tau= 0 : 0.005 : (h-0.002)
        tau_index = tau_index +1;
        lambda = eig(A);
        F_01 = inv_A* expm(A*h) * B;
        Av = inv_A * V;
        v_invB = inv(V) * B;
        
        F_11 = [Av(1,1)*v_invB(1,1); Av(2,1)*v_invB(1,1)];
        F_21 = [Av(1,2)*v_invB(2,1); Av(2,2)*v_invB(2,1)];
        F0 = [expm(A*h) F_01; 0 0 0];
        F1 = [zeros(2,2) F_11; zeros(1,3)];
        F2 = [zeros(2,2) F_21; zeros(1,3)];

        % the smaller the tau the bigger the alpha
        alpha1 = exp(lambda(1)*(h-tau));
        alpha2 = exp(lambda(2) * (h-tau));
        % Make them into vectors
        alpha2_v(h_index, tau_index) = alpha2;
        alpha1_v(h_index, tau_index) = alpha1;
        
        % find the minima and maxima of alphas
        [alpha1_max, I1_max] = max(alpha1_v(h_index,1:tau_index));
        [alpha1_min, I1_min] = min(alpha1_v(h_index,1:tau_index));
        [alpha2_max, I2_max] = max(alpha2_v(h_index,1:tau_index));
        [alpha2_min, I2_min] = min(alpha2_v(h_index,1:tau_index));

        F(:,:,h_index, tau_index) = F0 - alpha1 * F1 - alpha2 * F2;
        F23(h_index, tau_index) = F(2,3,h_index,tau_index);
        F13(h_index, tau_index) = F(1,3,h_index, tau_index);
        
        G0 = [-inv_A*B; 1];
        G1 = [F_11; 0];
        G2 = [F_21; 0];
        G(:,:,h_index, tau_index) = G0 + alpha1 * G1 + alpha2 * G2;
        K2(:,:,h_index, tau_index) = K_cl2;

    end
end

for x =5:5:50
    [alpha1_max, I11_max] = max(alpha1_v(x,(1:ceil(1*x/6))));
    [alpha1_min, I11_min] = min(alpha1_v(x,(1:ceil(1*x/6))));
    [alpha2_max, I12_max] = max(alpha2_v(x,(1:ceil(1*x/6))));
    [alpha2_min, I12_min] = min(alpha2_v(x,(1:ceil(1*x/6))));

    [alpha10_max, I10_max] = max(alpha1_v(x,(ceil(x/6):ceil(2*x/6))));
    [alpha10_min, I10_min] = min(alpha1_v(x,(ceil(x/6):ceil(2*x/6))));
    [alpha20_max, I20_max] = max(alpha2_v(x,(ceil(x/6):ceil(2*x/6))));
    [alpha20_min, I20_min] = min(alpha2_v(x,(ceil(x/6):ceil(2*x/6))));
        I10_max = I10_max + floor(x/6);
        I10_min = I10_min + floor(x/6);
        I20_max = I20_max + floor(x/6);
        I20_min = I20_min + floor(x/6);

    [alpha1_max, I21_max] = max(alpha1_v(x,(ceil(2*x/6):ceil(3*x/6))));
    [alpha1_min, I21_min] = min(alpha1_v(x,(ceil(2*x/6):ceil(3*x/6))));
    [alpha2_max, I22_max] = max(alpha2_v(x,(ceil(2*x/6):ceil(3*x/6))));
    [alpha2_min, I22_min] = min(alpha2_v(x,(ceil(2*x/6):ceil(3*x/6))));
        I21_max = I21_max + ceil(2*x/6)-1;
        I21_min = I21_min + ceil(2*x/6)-1;
        I22_max = I22_max + ceil(2*x/6)-1;
        I22_min = I22_min + ceil(2*x/6)-1;

    [alpha1_max, I31_max] = max(alpha1_v(x,(ceil(3*x/6):ceil(4*x/6))));
    [alpha1_min, I31_min] = min(alpha1_v(x,(ceil(3*x/6):ceil(4*x/6))));
    [alpha2_max, I32_max] = max(alpha2_v(x,(ceil(3*x/6):ceil(4*x/6))));
    [alpha2_min, I32_min] = min(alpha2_v(x,(ceil(3*x/6):ceil(4*x/6))));
        I31_max = I31_max + ceil(3*x/6)-1;
        I31_min = I31_min + ceil(3*x/6)-1;
        I32_max = I32_max + ceil(3*x/6)-1;
        I32_min = I32_min + ceil(3*x/6)-1;

    [alpha1_max, I41_max] = max(alpha1_v(x,ceil(4*x/6):ceil(5*x/6)));
    [alpha1_min, I41_min] = min(alpha1_v(x,ceil(4*x/6):ceil(5*x/6)));
    [alpha2_max, I42_max] = max(alpha2_v(x,ceil(4*x/6):ceil(5*x/6)));
    [alpha2_min, I42_min] = min(alpha2_v(x,ceil(4*x/6):ceil(5*x/6)));
        I41_max = I41_max + ceil(4*x/6)-1;
        I41_min = I41_min + ceil(4*x/6)-1;
        I42_max = I42_max + ceil(4*x/6)-1; 
        I42_min = I42_min + ceil(4*x/6)-1;

    [alpha1_max, I51_max] = max(alpha1_v(x,ceil(5*x/6):x));
    [alpha1_min, I51_min] = min(alpha1_v(x,1:x));
    [alpha2_max, I52_max] = max(alpha2_v(x,ceil(5*x/6):x));
    [alpha2_min, I52_min] = min(alpha2_v(x,1:x));
        I51_max = I51_max + ceil(5*x/6)-1;
%         I51_min = I51_min + I51_max -1 
        I52_max = I52_max + ceil(5*x/6)-1;
%         I52_min = I52_min + ceil(5*x/6);

        %mini Square
    Xs = [F23(x,I11_min), F23(x, I11_min), F23(x, I11_max), F23(x, I11_max), F23(x,I21_min), F23(x, I21_min), F23(x, I21_max), F23(x, I21_max), F23(x,I31_min), F23(x, I31_min), F23(x, I31_max), F23(x, I31_max), F23(x,I41_min), F23(x, I41_min), F23(x, I41_max), F23(x, I41_max) ];
    Ys = [F13(x, I12_min), F13(x, I12_max), F13(x, I12_max),F13(x, I12_min),   F13(x, I22_min), F13(x, I22_max), F13(x, I22_max),F13(x, I22_min), F13(x, I32_min), F13(x, I32_max), F13(x, I32_max),F13(x, I32_min), F13(x, I42_min), F13(x, I42_max), F13(x, I42_max),F13(x, I42_min)];
    
    Xs1 = [F23(x, I11_max)];%, F23(x, I11_max)];
    Ys1 = [F13(x, I12_max)];%,F13(x, I12_min)];
    Xs11 = [F23(x, I11_max), F23(x, I11_max)];
    Ys11 = [F13(x, I12_max),F13(x, I12_min)];
    Xs22 = [F23(x,I21_min), F23(x, I21_min), F23(x, I21_max), F23(x, I21_max)];
    Ys22 = [F13(x, I22_min), F13(x, I22_max), F13(x, I22_max),F13(x, I22_min)];
    % new polytope
    Xs2 = [F23(x, I21_max)];
    Ys2 = [F13(x, I22_min)];
    Xs33 = [F23(x,I31_min), F23(x, I31_min), F23(x, I31_max), F23(x, I31_max)];
    Ys33 = [F13(x, I32_min), F13(x, I32_max), F13(x, I32_max),F13(x, I32_min)];
    Xs3 = [F23(x, I31_max)];
    Ys3 = [F13(x, I32_min)];
    Xs44 = [F23(x,I41_min), F23(x, I41_min), F23(x, I41_max), F23(x, I41_max)];
    Ys44 = [F13(x, I42_min), F13(x, I42_max), F13(x, I42_max),F13(x, I42_min)];
%     Xs4 = [F23(x, I41_max)];
%     Ys4 = [F13(x, I42_min)];
    Xs55 = [F23(x,I10_min), F23(x, I10_min), F23(x, I10_max), F23(x, I10_max)];
    Ys55 = [F13(x, I20_min), F13(x, I20_max), F13(x, I20_max),F13(x, I20_min)];
%     Xs5 = [F23(x, I10_max)];
%     Ys5 = [F13(x, I20_min)];
    Xs66 = [F23(x,I51_min), F23(x, I51_min), F23(x, I51_max), F23(x, I51_max)];
    Ys66 = [F13(x, I52_min), F13(x, I52_max), F13(x, I52_max),F13(x, I52_min)];
    Xs6 = [F23(x, I51_max), F23(x,I51_min)];
    Ys6 = [F13(x, I52_min), F13(x, I52_min)];

    X = [Xs1  Xs2 Xs3  Xs6];
    Y = [Ys1  Ys2 Ys3  Ys6];

    I1_mean = mean(alpha1_v(x,(1:1*x)));
    [alpha1_max, I1_max] = max(alpha1_v(x,(1:1*x)));
    [alpha1_min, I1_min] = min(alpha1_v(x,(1:1*x)));
    [alpha2_max, I2_max] = max(alpha2_v(x,(1:1*x)));
    [alpha2_min, I2_min] = min(alpha2_v(x,(1:1*x)));
    I2_mean = mean(alpha2_v(x,(1:1*x)));
    [~, idx] = min(abs(alpha1_v(x,(1:1*x)) - I1_mean));
    [~, idx2] = min(abs(alpha2_v(x,(1:1*x)) - I2_mean));

    % Triangle
    Xt = [F23(x,I1_min),  F23(x, idx), F23(x, I1_max) ];
    Yt = [F13(x, I2_min),  F13(x, I2_min),F13(x, I2_max)];

    % Small polytope
    Xp = [F23(x,I1_min), F23(x,idx), F23(x, I1_max)];
    Yp = [F13(x, I2_min), F13(x, I2_min),  F13(x, I2_max)];


    figure(x+2)
    plot(F23(x,1:1*x), F13(x,1:1*x))
    hold on
    % title('Small Polytope')
    fill(X,Y,'r','FaceAlpha',0.2)
    scatter(X,Y, 'filled')
    fill(Xt,Yt,'y','FaceAlpha',0.2)
    scatter(Xt,Yt, 'filled')
    %     fill(Xs11,Ys11,'y','FaceAlpha',0.2)
    % scatter(Xs11,Ys11, 'filled')
    %     fill(Xs22,Ys22,'y','FaceAlpha',0.2)
    % scatter(Xs22,Ys22, 'filled')
    %     fill(Xs33,Ys33,'y','FaceAlpha',0.2)
    % scatter(Xs33,Ys33, 'filled')
    %     fill(Xs44,Ys44,'y','FaceAlpha',0.2)
    % scatter(Xs44,Ys44, 'filled')
    %     fill(Xs55,Ys55,'y','FaceAlpha',0.2)
    % scatter(Xs55,Ys55, 'filled')
    %     fill(Xs66,Ys66,'y','FaceAlpha',0.2)
    % scatter(Xs66,Ys66, 'filled')
    xlabel('$F_{2,3}$', 'Interpreter', ' latex')
    ylabel('$F_{1,3}$', 'Interpreter', ' latex')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 3b
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve LMIs for stability anlaysis
[V,D] = eig(A);
inv_A = inv(A);
h_index = 0;
h_stable_LMI = [];
tau_stable_LMI = [];
tau_UNstable_LMI =[];
h_UNstable_LMI = [];
for h = 0.0005: 0.0005:0.5
%     h = 0.05
    h_index = h_index +1;
    tau_index = 0;
    % Compute only for h 
    lambda = eig(A);
    F_01 = inv_A* expm(A*h) * B;
    Av = inv_A * V;
    v_invB = inv(V) * B;
    F_11 = [Av(1,1)*v_invB(1,1); Av(2,1)*v_invB(1,1)];
    F_21 = [Av(1,2)*v_invB(2,1); Av(2,2)*v_invB(2,1)];
    F0 = [expm(A*h) F_01; 0 0 0];
    F1 = [zeros(2,2) F_11; zeros(1,3)];
    F2 = [zeros(2,2) F_21; zeros(1,3)];
    G0 = [-inv_A*B; 1];
    G1 = [F_11; 0];
    G2 = [F_21; 0];

    for tau= 0 : 0.0005 : (h-0.005)
        tau_index = tau_index +1; 
%         tau = 0.01;
        % define alpha vectors
        alpha1 = exp(lambda(1)*(h-tau));
        alpha2 = exp(lambda(2) * (h-tau));
        alpha2_v(h_index, tau_index) = alpha2;
        alpha1_v(h_index, tau_index) = alpha1;
        % find the minima and maxima of alphas
        [alpha1_max, I1_max] = max(alpha1_v(h_index,1:tau_index));
        [alpha1_min, I1_min] = min(alpha1_v(h_index,1:tau_index));
        [alpha2_max, I2_max] = max(alpha2_v(h_index,1:tau_index));
        [alpha2_min, I2_min] = min(alpha2_v(h_index,1:tau_index));
        
%         cl_compare = (F0 - alpha1 * F1 - alpha2 * F2)-(G0 + alpha1 * G1 + alpha2 * G2)*K_cl2;
%         F_compare =    (F0 - alpha1 * F1 - alpha2 * F2)
%         G_compare =  (G0 + alpha1 * G1 + alpha2 * G2)
         % Compute matrix at the vertices: F
         F1_min = F0 - alpha1_min * F1 - alpha2_min * F2;
         F2_min = F0 - alpha1_max * F1 - alpha2_max * F2;
         F3_min = F0 - alpha1_max * F1 - alpha2_min * F2;
         % G
         G1_min = G0 + alpha1_min * G1 + alpha2_min * G2;
         G2_min = G0 + alpha1_max * G1 + alpha2_max * G2;
         G3_min = G0 + alpha1_max * G1 + alpha2_min * G2;


         % Closed loop matrices: H
         H_h11 = F1_min - G1_min * K_cl2;
         H_h21 = F2_min - G1_min * K_cl2;
         H_h31 = F3_min - G1_min * K_cl2;

         % LMI
            gamma = 0.2;
            P = sdpvar(3,3);
            Q = eye(3);
            cons = [P >= eye(3)*0.0002, 
                H_h11' * P * H_h11 - P <= - Q,
                H_h21' * P * H_h21 - P <= - Q,
                H_h31' * P * H_h31 - P <= - Q,
                ];
    
            obj = 0;
            S = sdpsettings('solver', 'sedumi');
            result = optimize(cons,obj, S);
%             disp(result.info)
%             value(P)
    
            if (result.problem == 0) 
                h_stable_LMI(end+1) = h;
                tau_stable_LMI(end+1) = tau;
            else
                h_UNstable_LMI(end+1) = h;
                tau_UNstable_LMI(end+1) = tau;
    
            end
     end
end


%%
% figure(32)
% scatter(F1_min(2,3), F1_min(1,3))
% hold on
% scatter(F2_min(2,3), F2_min(1,3))
% scatter(F3_min(2,3), F3_min(1,3))


figure(31)
scatter(h_stable_LMI, tau_stable_LMI, 'green' ) 
hold on
scatter(h_UNstable_LMI, tau_UNstable_LMI,'red' ) 
xlabel('Sampling Interval (h)');
ylabel('System Delay (\tau)');
title('Stable delays and sampling time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Question 3c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Solve LMIs for stability anlaysis
[V,D] = eig(A);
inv_A = inv(A);
h_index = 0;
h_stable_LMI = [];
tau_stable_LMI = [];
tau_UNstable_LMI =[];
h_UNstable_LMI = [];
x =0;
for h = 0.0005: 0.0005:0.1
%     h = 0.05
    h_index = h_index +1;
    tau_index = 0;
    % Compute only for h 
    lambda = eig(A);
    F_01 = inv_A* expm(A*h) * B;
    Av = inv_A * V;
    v_invB = inv(V) * B;
    F_11 = [Av(1,1)*v_invB(1,1); Av(2,1)*v_invB(1,1)];
    F_21 = [Av(1,2)*v_invB(2,1); Av(2,2)*v_invB(2,1)];
    F0 = [expm(A*h) F_01; 0 0 0];
    F1 = [zeros(2,2) F_11; zeros(1,3)];
    F2 = [zeros(2,2) F_21; zeros(1,3)];
    G0 = [-inv_A*B; 1];
    G1 = [F_11; 0];
    G2 = [F_21; 0];

    % tau = 0.01;
    for tau= 0 : 0% 0.0005 : (h-0.0005)
        tau_index = tau_index +1;  
        % defin alpha vectors
        alpha1 = exp(lambda(1)*(h-tau));
        alpha2 = exp(lambda(2) * (h-tau));
        alpha2_v(h_index, tau_index) = alpha2;
        alpha1_v(h_index, tau_index) = alpha1;

        % find the minima and maxima of alphas
        [alpha1_max, I1_max] = max(alpha1_v(h_index,1:tau_index));
        [alpha1_min, I1_min] = min(alpha1_v(h_index,1:tau_index));
        [alpha2_max, I2_max] = max(alpha2_v(h_index,1:tau_index));
        [alpha2_min, I2_min] = min(alpha2_v(h_index,1:tau_index));
        
        I1_mean = mean(alpha1_v(h_index,1:tau_index));
        I2_mean = mean(alpha2_v(h_index,1:tau_index));
        [~, idx] = min(abs(alpha1_v(h_index,1:tau_index) - I1_mean));
        [~, idx2] = min(abs(alpha2_v(h_index,1:tau_index) - I2_mean));
        alpha1_mean = alpha1_v(h_index, idx);
        alpha2_mean = alpha2_v(h_index, idx2);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     x = h_index;
    [alpha11_max, I11_max] = max(alpha1_v(x,(1:ceil(1*x/6))));
    [alpha12_min, I11_min] = min(alpha1_v(x,(1:ceil(1*x/6))));
    [alpha12_max, I12_max] = max(alpha2_v(x,(1:ceil(1*x/6))));
    [alpha12_min, I12_min] = min(alpha2_v(x,(1:ceil(1*x/6))));

    [alpha21_max, I21_max] = max(alpha1_v(x,(ceil(2*x/6):ceil(3*x/6))));
    [alpha21_min, I21_min] = min(alpha1_v(x,(ceil(2*x/6):ceil(3*x/6))));
    [alpha22_max, I22_max] = max(alpha2_v(x,(ceil(2*x/6):ceil(3*x/6))));
    [alpha22_min, I22_min] = min(alpha2_v(x,(ceil(2*x/6):ceil(3*x/6))));

    [alpha31_max, I31_max] = max(alpha1_v(x,(ceil(3*x/6):ceil(4*x/6))));
    [alpha31_min, I31_min] = min(alpha1_v(x,(ceil(3*x/6):ceil(4*x/6))));
    [alpha32_max, I32_max] = max(alpha2_v(x,(ceil(3*x/6):ceil(4*x/6))));
    [alpha32_min, I32_min] = min(alpha2_v(x,(ceil(3*x/6):ceil(4*x/6))));

    [alpha51_max, I51_max] = max(alpha1_v(x,ceil(5*x/6):x));
    [alpha51_min, I51_min] = min(alpha1_v(x,1:x));
    [alpha52_max, I52_max] = max(alpha2_v(x,ceil(5*x/6):x));
    [alpha52_min, I52_min] = min(alpha2_v(x,1:x));
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Small polytope
%         Xp1 = [F23(x,I1_min), F23(x,idx), F23(x, I1_max), F23(x, I1_max)];
%         Yp1 = [F13(x, I2_min), F13(x, I2_min), F13(x, idx2), F13(x, I2_max)];
            
        % Compute matrix at the vertices: F
        F1_min = F0 - alpha11_max * F1 - alpha12_max * F2;
        F2_min = F0 - alpha21_max * F1 - alpha22_min * F2;
        F3_min = F0 - alpha31_max * F1 - alpha32_min * F2;
        F4_min = F0 - alpha51_max * F1 - alpha52_min * F2;
        F5_min = F0 - alpha51_min * F1 - alpha52_min * F2;

        % G
        G1_min = G0 - alpha11_max * G1 - alpha12_max * G2;
        G2_min = G0 - alpha21_max * G1 - alpha22_min * G2;
        G3_min = G0 - alpha31_max * G1 - alpha32_min * G2;
        G4_min = G0 - alpha51_max * G1 - alpha52_min * G2;
        G5_min = G0 - alpha51_min * G1 - alpha52_min * G2;
        

         % Closed loop matrices: H
        H_h1 = F1_min - G1_min * K_cl2;
        H_h2 = F2_min - G2_min * K_cl2;
        H_h3 = F3_min - G3_min * K_cl2;
        H_h4 = F4_min - G4_min * K_cl2;
        H_h5 = F5_min - G5_min * K_cl2;

         % LMI
            P = sdpvar(3,3);
            gamma = 0.2;
            Q = eye(3);
            cons = [P >= eye(3)*0.002, 
%                 H_h1' * P * H_h1 - P <= - Q,
%                 H_h2' * P * H_h2 - P <= - Q,
%                 H_h3' * P * H_h3 - P <= - Q,
                H_h4' * P * H_h4 - P <= - Q,
%                 H_h5' * P * H_h5 - P <= - Q,
                ];  
            obj = 0;
            S = sdpsettings('solver', 'sedumi');
            result = optimize(cons,obj, S);
            if (result.problem == 0) 
                h_stable_LMI(end+1) = h;
                tau_stable_LMI(end+1) = tau;
            else
                h_UNstable_LMI(end+1) = h;
                tau_UNstable_LMI(end+1) = tau;
    
            end
    end
end
%%

figure(31)
scatter(h_stable_LMI, tau_stable_LMI, 'green' ) 
hold on
scatter(h_UNstable_LMI, tau_UNstable_LMI,'red' ) 
xlabel('Sampling Interval (h)');
ylabel('System Delay (\tau)');
title('Stable delays and sampling time')