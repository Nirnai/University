function [ tau ] = Tau_OperationalSpace( u )
%TAU Summary of this function goes here
%   Detailed explanation goes here
% Times are defined in TrajGen.m file


% Robot Parameters
% m2=u(39); 
% m3=u(40); 
% g=u(41); 
% L1=u(42); 
% L2=u(43); 
% L3=u(44); 
% L4=u(45); 
% L5=u(46); 
% L6=u(47); 
% L7=u(48); 
% L8=u(49); 
% L9=u(50); 
% L10=u(51); 
% gx=u(52); 
% gy=u(53); 
% gz=u(54);
% I133=u(60); 
% I211=u(61); 
% I212=u(62); 
% I213=u(63); 
% I222=u(64); 
% I223=u(65); 
% I233=u(66); 
% I311=u(67); 
% I312=u(68);
% I313=u(69); 
% I322=u(70); 
% I323=u(71); 
% I333=u(72);
% 
% 
% 
% 
% t=u(1);
% 
% 
% 
% 
% 
% %% State
% q1=u(23);
% q2=u(24);
% q3=u(25);
% 
% Q = [q1;q2;q3];
% 
% qp1=u(26);
% qp2=u(27);
% qp3=u(28);
% 
% Qp = [qp1;qp2;qp3];
% 
% 
% 
% 
% 
% %% Jacobian
% 
% J =[ ...
%   cos(q1)*(L2 + L4) + L3*cos(q2)*sin(q1) + L5*cos(q2)*cos(q3)*sin(q1) - L5*sin(q1)*sin(q2)*sin(q3),                                                                                                                                                                                  cos(q1)*(L5*sin(q2 + q3) + L3*sin(q2)),                                                                                                                                                       L5*sin(q2 + q3)*cos(q1);
%   sin(q1)*(L2 + L4) - L3*cos(q1)*cos(q2) - L5*cos(q1)*cos(q2)*cos(q3) + L5*cos(q1)*sin(q2)*sin(q3),                                                                                                                                                                                  sin(q1)*(L5*sin(q2 + q3) + L3*sin(q2)),                                                                                                                                                       L5*sin(q2 + q3)*sin(q1);
%                                                                                                  0, cos(q1)*(sin(q1)*(L2 + L4) - L3*cos(q1)*cos(q2) - L5*cos(q1)*cos(q2)*cos(q3) + L5*cos(q1)*sin(q2)*sin(q3)) - sin(q1)*(cos(q1)*(L2 + L4) + L3*cos(q2)*sin(q1) + L5*cos(q2)*cos(q3)*sin(q1) - L5*sin(q1)*sin(q2)*sin(q3)), cos(q1)*(sin(q1)*(L2 + L4) - L5*cos(q1)*cos(q2)*cos(q3) + L5*cos(q1)*sin(q2)*sin(q3)) - sin(q1)*(cos(q1)*(L2 + L4) + L5*cos(q2)*cos(q3)*sin(q1) - L5*sin(q1)*sin(q2)*sin(q3));
%                                                                                                  0,                                                                                                                                                                                                                 sin(q1),                                                                                                                                                                       sin(q1);
%                                                                                                  0,                                                                                                                                                                                                                -cos(q1),                                                                                                                                                                      -cos(q1);
%                                                                                                  1,                                                                                                                                                                                                                       0,                                                                                                                                                                             0];
%  
%  
% 
% Jp = [... 
%   - qp1*(sin(q1)*(L2 + L4) - L3*cos(q1)*cos(q2) - L5*cos(q1)*cos(q2)*cos(q3) + L5*cos(q1)*sin(q2)*sin(q3)) - qp2*(L3*sin(q1)*sin(q2) + L5*cos(q2)*sin(q1)*sin(q3) + L5*cos(q3)*sin(q1)*sin(q2)) - qp3*(L5*cos(q2)*sin(q1)*sin(q3) + L5*cos(q3)*sin(q1)*sin(q2)),                                                                                                                                                                                                               qp2*cos(q1)*(L5*cos(q2 + q3) + L3*cos(q2)) - qp1*sin(q1)*(L5*sin(q2 + q3) + L3*sin(q2)) + L5*qp3*cos(q2 + q3)*cos(q1),                                                                                                                                                                                                   L5*qp2*cos(q2 + q3)*cos(q1) + L5*qp3*cos(q2 + q3)*cos(q1) - L5*qp1*sin(q2 + q3)*sin(q1);
%     qp2*(L3*cos(q1)*sin(q2) + L5*cos(q1)*cos(q2)*sin(q3) + L5*cos(q1)*cos(q3)*sin(q2)) + qp1*(cos(q1)*(L2 + L4) + L3*cos(q2)*sin(q1) + L5*cos(q2)*cos(q3)*sin(q1) - L5*sin(q1)*sin(q2)*sin(q3)) + qp3*(L5*cos(q1)*cos(q2)*sin(q3) + L5*cos(q1)*cos(q3)*sin(q2)),                                                                                                                                                                                                               qp2*sin(q1)*(L5*cos(q2 + q3) + L3*cos(q2)) + qp1*cos(q1)*(L5*sin(q2 + q3) + L3*sin(q2)) + L5*qp3*cos(q2 + q3)*sin(q1),                                                                                                                                                                                                   L5*qp1*sin(q2 + q3)*cos(q1) + L5*qp2*cos(q2 + q3)*sin(q1) + L5*qp3*cos(q2 + q3)*sin(q1);
%                                                                                                                                                                                                                                                               0, qp2*(cos(q1)*(L3*cos(q1)*sin(q2) + L5*cos(q1)*cos(q2)*sin(q3) + L5*cos(q1)*cos(q3)*sin(q2)) + sin(q1)*(L3*sin(q1)*sin(q2) + L5*cos(q2)*sin(q1)*sin(q3) + L5*cos(q3)*sin(q1)*sin(q2))) + qp3*(cos(q1)*(L5*cos(q1)*cos(q2)*sin(q3) + L5*cos(q1)*cos(q3)*sin(q2)) + sin(q1)*(L5*cos(q2)*sin(q1)*sin(q3) + L5*cos(q3)*sin(q1)*sin(q2))), qp2*(cos(q1)*(L5*cos(q1)*cos(q2)*sin(q3) + L5*cos(q1)*cos(q3)*sin(q2)) + sin(q1)*(L5*cos(q2)*sin(q1)*sin(q3) + L5*cos(q3)*sin(q1)*sin(q2))) + qp3*(cos(q1)*(L5*cos(q1)*cos(q2)*sin(q3) + L5*cos(q1)*cos(q3)*sin(q2)) + sin(q1)*(L5*cos(q2)*sin(q1)*sin(q3) + L5*cos(q3)*sin(q1)*sin(q2)));
%                                                                                                                                                                                                                                                               0,                                                                                                                                                                                                                                                                                                                         qp1*cos(q1),                                                                                                                                                                                                                                                                               qp1*cos(q1);
%                                                                                                                                                                                                                                                               0,                                                                                                                                                                                                                                                                                                                         qp1*sin(q1),                                                                                                                                                                                                                                                                               qp1*sin(q1);
%                                                                                                                                                                                                                                                               0,                                                                                                                                                                                                                                                                                                                                   0,                                                                                                                                                                                                                                                                                         0];
%                                                                                       
%                                                                                                   
%                                                                                                   
% X = [...
%  sin(q1)*(L2 + L4) - L3*cos(q1)*cos(q2) - L5*cos(q1)*cos(q2)*cos(q3) + L5*cos(q1)*sin(q2)*sin(q3);
%  L5*sin(q1)*sin(q2)*sin(q3) - L3*cos(q2)*sin(q1) - L5*cos(q2)*cos(q3)*sin(q1) - cos(q1)*(L2 + L4);
%                                                                 L1 - L5*sin(q2 + q3) - L3*sin(q2)];
%                                                             
% 
%                                                             
% w = det(sqrt(J'*J));
% %% Operational State
% Xp = J(1:3,:) * Qp;                                                                                                     
%                                                                                                       
% %% Desired State
% a = u(2);
% b = u(3);
% c = u(4);
% 
% Q_Desired = u(5:13);
% X_Desired = u(14:22);    
% 
% Qd = Q_Desired(1:3);
% Qdp = Q_Desired(4:6);
% Qdpp = Q_Desired(7:9);
% 
% Xd = X_Desired(1:3);
% Xdp = X_Desired(4:6);
% Xdpp = X_Desired(7:9);
% 
% %% Controller Parameters
% Kp=diag([u(29);u(30);u(31)]);
% Kd=diag([u(32);u(33);u(34)]);
% Ki=diag([u(35);u(36);u(37)]);
% 
% 
% %% Errors
% if(a)
%     DeltaQ = Q - Qd;
%     DeltaQp = Qp - Qdp;
% 
%     Qrp = Qdp - Kp * DeltaQ;
%     Qrpp = Qdpp - Kp * DeltaQp;
%     
% elseif(b)
%     % Error in Position
%     DeltaX = X - Xd;
%     % Error in Velocity
%     DeltaXp = Xp - Xdp;
% 
%     % PD-Like Dynamics
%     Xrp = Xdp - Kp * DeltaX;
%     Xrpp = Xdpp - Kp * DeltaXp;
%     Qrp = J(1:3,:)\Xrp;
%     Qrpp = J(1:3,:)\(Xrpp - Jp(1:3,:)*Qp);
% %     % Integrated Error in Position
% %     persistent DeltaXTotal previousTime
% %     if(t == 0)
% %         DeltaXTotal = DeltaX;
% %         previousTime = t;
% %     else
% %         deltaT = t - previousTime;
% %         DeltaXTotal = DeltaXTotal + (DeltaX * deltaT);
% %         previousTime = t;
% %     end
% % 
% elseif(c)
% %     % Error in Position
% %     DeltaX = X - Xd;
% %     % Error in Velocity
% %     DeltaXp = Xp - Xdp;
% % 
% %     % PD-Like Dynamics
% %     Xrp = Xdp - Kp * DeltaX;
% %     Xrpp = Xdpp - Kp * DeltaXp;
% %     Qrp = J(1:3,:)\Xrp;
% %     Qrpp = J(1:3,:)\(Xrpp - Jp(1:3,:)*Qp);
%     Qrp = zeros(3,1);
%     Qrpp = zeros(3,1);
%     
% end
% 
% 
% %% Desired Operational Velocity Dynamics
% 
% % Break 
% % Xrp = [0;0;0];
% % Xrpp = [0;0;0];
% 
% % PD-Like Dynamics
% % Xrp = Xdp - Kp * DeltaX;
% % Xrpp = Xdpp - Kp * DeltaXp;
% 
% % PD-Like Dynamics
% % Xrp = Xdp - Kp * DeltaX - Ki * DeltaXTotal;
% % Xrpp = Xdpp - Kp * DeltaXp - Ki * DeltaX;
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% Desired Join Velocity Dynamics
% % Qrp = J(1:3,:)\Xrp;
% % Qrpp = J(1:3,:)\(Xrpp - Jp(1:3,:)*Qp);
% 
% qrp1 = Qrp(1);
% qrp2 = Qrp(2);
% qrp3 = Qrp(3);
% 
% qrpp1 = Qrpp(1);
% qrpp2 = Qrpp(2);
% qrpp3 = Qrpp(3);
% 
% %% Error from desired dynamics
% Sq = Qp - Qrp;
% 
% 
% %% Robot Desired Dynamics
% 
% % Robot Parameter Vector
% Theta = [...   
%      L3^2*m3;
%     L3*L9*m3;
%    L3*L10*m3;
%   L3*g*gx*m3;
%   L3*g*gy*m3;
%      L7^2*m2;
%     L7*L8*m2;
%   L7*g*gx*m2;
%   L7*g*gy*m2;
%      L8^2*m2;
%   L8*g*gx*m2;
%   L8*g*gy*m2;
%      L9^2*m3;
%    L9*L10*m3;
%   L9*g*gx*m3;
%   L9*g*gy*m3;
%     L10^2*m3;
%  L10*g*gx*m3;
%  L10*g*gy*m3;
%         I133;
%         I211;
%         I212;
%         I213;
%         I222;
%         I223;
%         I311;
%         I312;
%         I313;
%         I322;
%         I323;
%   L3*g*gz*m3;
%   L8*g*gz*m2;
%         I233;
%  L10*g*gz*m3;
%         I333];
%      
% % Robot Regressor
% Yr = [...
%   qrpp1*(cos(2*q2)/2 + 1/2) - (qp1*qrp2*sin(2*q2))/2 - (qp2*qrp1*sin(2*q2))/2, qrpp2*sin(q2) + qp2*qrp2*cos(q2), qrpp1*(cos(2*q2 + q3) + cos(q3)) - qrp1*((qp3*(sin(2*q2 + q3) + sin(q3)))/2 + qp2*sin(2*q2 + q3)) - qp1*qrp2*sin(2*q2 + q3) - (qp1*qrp3*(sin(2*q2 + q3) + sin(q3)))/2, cos(q2)*sin(q1), -cos(q1)*cos(q2), qrpp1, qrpp2*sin(q2) + qp2*qrp2*cos(q2), cos(q1), sin(q1), qrpp1*(cos(2*q2)/2 + 1/2) - (qp1*qrp2*sin(2*q2))/2 - (qp2*qrp1*sin(2*q2))/2, cos(q2)*sin(q1), -cos(q1)*cos(q2), qrpp1, qrp2*(qp2*cos(q2 + q3) + qp3*cos(q2 + q3)) + qrp3*(qp2*cos(q2 + q3) + qp3*cos(q2 + q3)) + qrpp2*sin(q2 + q3) + qrpp3*sin(q2 + q3), cos(q1), sin(q1), qrpp1*(cos(2*q2 + 2*q3)/2 + 1/2) - qrp1*((qp2*sin(2*q2 + 2*q3))/2 + (qp3*sin(2*q2 + 2*q3))/2) - (qp1*qrp2*sin(2*q2 + 2*q3))/2 - (qp1*qrp3*sin(2*q2 + 2*q3))/2, cos(q2)*cos(q3)*sin(q1) - sin(q1)*sin(q2)*sin(q3), cos(q1)*sin(q2)*sin(q3) - cos(q1)*cos(q2)*cos(q3), qrpp1, (qp1*qrp2*sin(2*q2))/2 - qrpp1*(cos(2*q2)/2 - 1/2) + (qp2*qrp1*sin(2*q2))/2, qrpp1*sin(2*q2) + qp1*qrp2*cos(2*q2) + qp2*qrp1*cos(2*q2), - qrpp2*sin(q2) - qp2*qrp2*cos(q2), qrpp1*(cos(2*q2)/2 + 1/2) - (qp1*qrp2*sin(2*q2))/2 - (qp2*qrp1*sin(2*q2))/2, qp2*qrp2*sin(q2) - qrpp2*cos(q2), qrp1*((qp2*sin(2*q2 + 2*q3))/2 + (qp3*sin(2*q2 + 2*q3))/2) - qrpp1*(cos(2*q2 + 2*q3)/2 - 1/2) + (qp1*qrp2*sin(2*q2 + 2*q3))/2 + (qp1*qrp3*sin(2*q2 + 2*q3))/2, qrpp1*sin(2*q2 + 2*q3) + qrp1*(qp2*cos(2*q2 + 2*q3) + qp3*cos(2*q2 + 2*q3)) + qp1*qrp2*cos(2*q2 + 2*q3) + qp1*qrp3*cos(2*q2 + 2*q3), - qrp2*(qp2*cos(q2 + q3) + qp3*cos(q2 + q3)) - qrp3*(qp2*cos(q2 + q3) + qp3*cos(q2 + q3)) - qrpp2*sin(q2 + q3) - qrpp3*sin(q2 + q3), qrpp1*(cos(2*q2 + 2*q3)/2 + 1/2) - qrp1*((qp2*sin(2*q2 + 2*q3))/2 + (qp3*sin(2*q2 + 2*q3))/2) - (qp1*qrp2*sin(2*q2 + 2*q3))/2 - (qp1*qrp3*sin(2*q2 + 2*q3))/2, qrp2*(qp2*sin(q2 + q3) + qp3*sin(q2 + q3)) + qrp3*(qp2*sin(q2 + q3) + qp3*sin(q2 + q3)) - qrpp2*cos(q2 + q3) - qrpp3*cos(q2 + q3),        0,        0,     0,             0,             0;
%                                                qrpp2 + (qp1*qrp1*sin(2*q2))/2,                    qrpp1*sin(q2),                                                               2*qrpp2*cos(q3) + qrpp3*cos(q3) - qp3*qrp2*sin(q3) + qp1*qrp1*sin(2*q2 + q3) - qrp3*sin(q3)*(qp2 + qp3), cos(q1)*sin(q2),  sin(q1)*sin(q2),     0,                    qrpp1*sin(q2),       0,       0,                                              qrpp2 + (qp1*qrp1*sin(2*q2))/2, cos(q1)*sin(q2),  sin(q1)*sin(q2),     0,                                                                                                                qrpp1*sin(q2 + q3),       0,       0,                                                                                                                 qrpp2 + qrpp3 + (qp1*qrp1*sin(2*q2 + 2*q3))/2, cos(q1)*cos(q2)*sin(q3) + cos(q1)*cos(q3)*sin(q2), cos(q2)*sin(q1)*sin(q3) + cos(q3)*sin(q1)*sin(q2),     0,                                                     -(qp1*qrp1*sin(2*q2))/2,                                       -qp1*qrp1*cos(2*q2),                     -qrpp1*sin(q2),                                                      (qp1*qrp1*sin(2*q2))/2,                   -qrpp1*cos(q2),                                                                                                                                -(qp1*qrp1*sin(2*q2 + 2*q3))/2,                                                                                                          -qp1*qrp1*cos(2*q2 + 2*q3),                                                                                                                 -qrpp1*sin(q2 + q3),                                                                                                                                 (qp1*qrp1*sin(2*q2 + 2*q3))/2,                                                                                                               -qrpp1*cos(q2 + q3), -cos(q2), -cos(q2), qrpp2, -cos(q2 + q3), qrpp2 + qrpp3;
%                                                                             0,                                0,                                                                                            qrpp2*cos(q3) + qp2*qrp2*sin(q3) + (qp1*qrp1*(sin(2*q2 + q3) + sin(q3)))/2,               0,                0,     0,                                0,       0,       0,                                                                           0,               0,                0,     0,                                                                                                                qrpp1*sin(q2 + q3),       0,       0,                                                                                                                 qrpp2 + qrpp3 + (qp1*qrp1*sin(2*q2 + 2*q3))/2, cos(q1)*cos(q2)*sin(q3) + cos(q1)*cos(q3)*sin(q2), cos(q2)*sin(q1)*sin(q3) + cos(q3)*sin(q1)*sin(q2),     0,                                                                           0,                                                         0,                                  0,                                                                           0,                                0,                                                                                                                                -(qp1*qrp1*sin(2*q2 + 2*q3))/2,                                                                                                          -qp1*qrp1*cos(2*q2 + 2*q3),                                                                                                                 -qrpp1*sin(q2 + q3),                                                                                                                                 (qp1*qrp1*sin(2*q2 + 2*q3))/2,                                                                                                               -qrpp1*cos(q2 + q3),        0,        0,     0, -cos(q2 + q3), qrpp2 + qrpp3];
%                                                                   


%% Controller 

% PD Values: 
%       Kp = [80,80,10];
%       Kd = [3,3,3];

% PID Values:
%       Kd = [1;1;1];
%       Kp = [50;50;50];
%       Ki = [600;600;600];

% tau = -Kd*Sq + Yr*Theta;
% tau = Yr*Theta;
tau = [0;0;0];
end

