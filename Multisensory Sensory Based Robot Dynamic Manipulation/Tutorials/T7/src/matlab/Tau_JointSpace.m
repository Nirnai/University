function [ tau ] = Tau_JointSpace( u )
%TAU Summary of this function goes here
%   Detailed explanation goes here
% Times are defined in TrajGen.m file

t=u(1);

qd1 = u(2);
qd2 = u(3);
qd3 = u(4);

Qd = [qd1;qd2;qd3];

qdp1 = u(5);
qdp2 = u(6);
qdp3 = u(7);

Qdp = [qdp1;qdp2;qdp3];

qdpp1 = u(8);
qdpp2 = u(9);
qdpp3 = u(10);

Qdpp = [qdpp1;qdpp2;qdpp3];


%% State
q1=u(11);
q2=u(12);
q3=u(13);

Q = [q1;q2;q3];

qp1=u(14);
qp2=u(15);
qp3=u(16);

Qp = [qp1;qp2;qp3];

%% Controller Parameters
Kp=diag([u(17);u(18);u(19)]);
Kd=diag([u(20);u(21);u(22)]);
Ki=diag([u(23);u(24);u(25)]);


%% Errors
% Error in Position
DeltaQ = Q - Qd;
% 
% % Integrated Error in Position
% persistent DeltaQTotal previousTime
% if(t == 0)
%     DeltaQTotal = DeltaQ;
%     previousTime = t;
% else
%     deltaT = t - previousTime;
%     DeltaQTotal = DeltaQTotal + (DeltaQ * deltaT);
%     previousTime = t;
% end
% 
% Error in Velocity
DeltaQp = Qp - Qdp;



%% Desired Velocity Dynamics

% Break
% Qrp = [0;0;0];
% Qrpp = [0;0;0];

% PD-Like Dynamics
% PD Parameters: 
% Kp = [20;30;20]
% Kd = [5;5;5]
Qrp = Qdp - Kp * DeltaQ;
Qrpp = Qdpp - Kp * DeltaQp;



% PID-Like Dynamics
% Qrp = Qdp - Kp * DeltaQ - Ki * DeltaQTotal;
% Qrpp = Qdpp - Kp * DeltaQp - Ki * DeltaQ;





qrp1 = Qrp(1);
qrp2 = Qrp(2);
qrp3 = Qrp(3);

qrpp1 = Qrpp(1);
qrpp2 = Qrpp(2);
qrpp3 = Qrpp(3);

%% Error from desired dynamics
Sq = Qp - Qrp;


%% Robot Desired Dynamics
% Robot Parameters
m2=u(27); m3=u(28); g=u(29); L3=u(32); L7=u(36); L8=u(37); L9=u(38); L10=u(39); gx=u(40); gy=u(41); gz=u(42);
I133=u(48); I211=u(49); I212=u(50); I213=u(51); I222=u(52); I223=u(53); I233=u(54); I311=u(55); I312=u(56);
I313=u(57); I322=u(58); I323=u(59); I333=u(60);

% Robot Parameter Vector
Theta = [...   
     L3^2*m3;
    L3*L9*m3;
   L3*L10*m3;
  L3*g*gx*m3;
  L3*g*gy*m3;
     L7^2*m2;
    L7*L8*m2;
  L7*g*gx*m2;
  L7*g*gy*m2;
     L8^2*m2;
  L8*g*gx*m2;
  L8*g*gy*m2;
     L9^2*m3;
   L9*L10*m3;
  L9*g*gx*m3;
  L9*g*gy*m3;
    L10^2*m3;
 L10*g*gx*m3;
 L10*g*gy*m3;
        I133;
        I211;
        I212;
        I213;
        I222;
        I223;
        I311;
        I312;
        I313;
        I322;
        I323;
  L3*g*gz*m3;
  L8*g*gz*m2;
        I233;
 L10*g*gz*m3;
        I333];
     
% Robot Regressor
Yr = [...
  qrpp1*(cos(2*q2)/2 + 1/2) - (qp1*qrp2*sin(2*q2))/2 - (qp2*qrp1*sin(2*q2))/2, qrpp2*sin(q2) + qp2*qrp2*cos(q2), qrpp1*(cos(2*q2 + q3) + cos(q3)) - qrp1*((qp3*(sin(2*q2 + q3) + sin(q3)))/2 + qp2*sin(2*q2 + q3)) - qp1*qrp2*sin(2*q2 + q3) - (qp1*qrp3*(sin(2*q2 + q3) + sin(q3)))/2, cos(q2)*sin(q1), -cos(q1)*cos(q2), qrpp1, qrpp2*sin(q2) + qp2*qrp2*cos(q2), cos(q1), sin(q1), qrpp1*(cos(2*q2)/2 + 1/2) - (qp1*qrp2*sin(2*q2))/2 - (qp2*qrp1*sin(2*q2))/2, cos(q2)*sin(q1), -cos(q1)*cos(q2), qrpp1, qrp2*(qp2*cos(q2 + q3) + qp3*cos(q2 + q3)) + qrp3*(qp2*cos(q2 + q3) + qp3*cos(q2 + q3)) + qrpp2*sin(q2 + q3) + qrpp3*sin(q2 + q3), cos(q1), sin(q1), qrpp1*(cos(2*q2 + 2*q3)/2 + 1/2) - qrp1*((qp2*sin(2*q2 + 2*q3))/2 + (qp3*sin(2*q2 + 2*q3))/2) - (qp1*qrp2*sin(2*q2 + 2*q3))/2 - (qp1*qrp3*sin(2*q2 + 2*q3))/2, cos(q2)*cos(q3)*sin(q1) - sin(q1)*sin(q2)*sin(q3), cos(q1)*sin(q2)*sin(q3) - cos(q1)*cos(q2)*cos(q3), qrpp1, (qp1*qrp2*sin(2*q2))/2 - qrpp1*(cos(2*q2)/2 - 1/2) + (qp2*qrp1*sin(2*q2))/2, qrpp1*sin(2*q2) + qp1*qrp2*cos(2*q2) + qp2*qrp1*cos(2*q2), - qrpp2*sin(q2) - qp2*qrp2*cos(q2), qrpp1*(cos(2*q2)/2 + 1/2) - (qp1*qrp2*sin(2*q2))/2 - (qp2*qrp1*sin(2*q2))/2, qp2*qrp2*sin(q2) - qrpp2*cos(q2), qrp1*((qp2*sin(2*q2 + 2*q3))/2 + (qp3*sin(2*q2 + 2*q3))/2) - qrpp1*(cos(2*q2 + 2*q3)/2 - 1/2) + (qp1*qrp2*sin(2*q2 + 2*q3))/2 + (qp1*qrp3*sin(2*q2 + 2*q3))/2, qrpp1*sin(2*q2 + 2*q3) + qrp1*(qp2*cos(2*q2 + 2*q3) + qp3*cos(2*q2 + 2*q3)) + qp1*qrp2*cos(2*q2 + 2*q3) + qp1*qrp3*cos(2*q2 + 2*q3), - qrp2*(qp2*cos(q2 + q3) + qp3*cos(q2 + q3)) - qrp3*(qp2*cos(q2 + q3) + qp3*cos(q2 + q3)) - qrpp2*sin(q2 + q3) - qrpp3*sin(q2 + q3), qrpp1*(cos(2*q2 + 2*q3)/2 + 1/2) - qrp1*((qp2*sin(2*q2 + 2*q3))/2 + (qp3*sin(2*q2 + 2*q3))/2) - (qp1*qrp2*sin(2*q2 + 2*q3))/2 - (qp1*qrp3*sin(2*q2 + 2*q3))/2, qrp2*(qp2*sin(q2 + q3) + qp3*sin(q2 + q3)) + qrp3*(qp2*sin(q2 + q3) + qp3*sin(q2 + q3)) - qrpp2*cos(q2 + q3) - qrpp3*cos(q2 + q3),        0,        0,     0,             0,             0;
                                               qrpp2 + (qp1*qrp1*sin(2*q2))/2,                    qrpp1*sin(q2),                                                               2*qrpp2*cos(q3) + qrpp3*cos(q3) - qp3*qrp2*sin(q3) + qp1*qrp1*sin(2*q2 + q3) - qrp3*sin(q3)*(qp2 + qp3), cos(q1)*sin(q2),  sin(q1)*sin(q2),     0,                    qrpp1*sin(q2),       0,       0,                                              qrpp2 + (qp1*qrp1*sin(2*q2))/2, cos(q1)*sin(q2),  sin(q1)*sin(q2),     0,                                                                                                                qrpp1*sin(q2 + q3),       0,       0,                                                                                                                 qrpp2 + qrpp3 + (qp1*qrp1*sin(2*q2 + 2*q3))/2, cos(q1)*cos(q2)*sin(q3) + cos(q1)*cos(q3)*sin(q2), cos(q2)*sin(q1)*sin(q3) + cos(q3)*sin(q1)*sin(q2),     0,                                                     -(qp1*qrp1*sin(2*q2))/2,                                       -qp1*qrp1*cos(2*q2),                     -qrpp1*sin(q2),                                                      (qp1*qrp1*sin(2*q2))/2,                   -qrpp1*cos(q2),                                                                                                                                -(qp1*qrp1*sin(2*q2 + 2*q3))/2,                                                                                                          -qp1*qrp1*cos(2*q2 + 2*q3),                                                                                                                 -qrpp1*sin(q2 + q3),                                                                                                                                 (qp1*qrp1*sin(2*q2 + 2*q3))/2,                                                                                                               -qrpp1*cos(q2 + q3), -cos(q2), -cos(q2), qrpp2, -cos(q2 + q3), qrpp2 + qrpp3;
                                                                            0,                                0,                                                                                            qrpp2*cos(q3) + qp2*qrp2*sin(q3) + (qp1*qrp1*(sin(2*q2 + q3) + sin(q3)))/2,               0,                0,     0,                                0,       0,       0,                                                                           0,               0,                0,     0,                                                                                                                qrpp1*sin(q2 + q3),       0,       0,                                                                                                                 qrpp2 + qrpp3 + (qp1*qrp1*sin(2*q2 + 2*q3))/2, cos(q1)*cos(q2)*sin(q3) + cos(q1)*cos(q3)*sin(q2), cos(q2)*sin(q1)*sin(q3) + cos(q3)*sin(q1)*sin(q2),     0,                                                                           0,                                                         0,                                  0,                                                                           0,                                0,                                                                                                                                -(qp1*qrp1*sin(2*q2 + 2*q3))/2,                                                                                                          -qp1*qrp1*cos(2*q2 + 2*q3),                                                                                                                 -qrpp1*sin(q2 + q3),                                                                                                                                 (qp1*qrp1*sin(2*q2 + 2*q3))/2,                                                                                                               -qrpp1*cos(q2 + q3),        0,        0,     0, -cos(q2 + q3), qrpp2 + qrpp3];
                                                                  


%% Controller 

% PD Values: 
%       Kp = [80,80,10];
%       Kd = [3,3,3];

% PID Values:

tau =-Kd*Sq + Yr*Theta;
% tau = Yr*Theta;

end

