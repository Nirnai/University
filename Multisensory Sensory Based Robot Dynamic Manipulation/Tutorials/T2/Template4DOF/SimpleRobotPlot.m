function [ Xout ] = SimpleRobotPlot( u )
%SIMPLEROBOTPLOT Summary of this function goes here

%%%%%FOLLOW THE INSTRUCTIONS GIVEN IN THE PDF FILE BEFORE STARTING THE
%%%%%IMPLEMENTATION

%Joint Position (NOT NEEDED IN THIS TUTORIAL)
qp1=u(1);
qp2=u(2);
qp3=u(3);
qp4=u(4);

%Joint Velocity
q1=u(5);
q2=u(6);
q3=u(7);
q4=u(8);

%Time
t=u(9);

%Kinematic Parameters
L1=u(10);
L2=u(11);
L3=u(12);
% ...
% TODO: add as many variable as needed, according to the robot parameters.

q1 = sym('q1');
q2 = sym('q2');
q3 = sym('q3');
q4 = sym('q4');

L4 = 0.08;
L5 = 0.08;
L6 = 0.115;
L7 = 0.082;
L8 = 0.108;
L9 = 0.123;
L10 = 0.167;
L11 = 0.028;
L13 = 0.148;
L14 = 0.062;

a1 = 0;
a2 = 0; 
a3 = 0.107; 
a4 = L13; 
a5 = 0; 
a6 = L7/2; 
a7 = L10/(sin(pi/4)); 
a8 = L13;

alpha1 = -pi/2;
alpha2 = pi/2;
alpha3 = 3/4 * pi;
alpha4 = 0;
alpha5 = 0;
alpha6 = 0;
alpha7 = 0;
alpha8 = 0;

theta1 = q1;
theta2 = pi/2;
theta3 = 0;
theta4 = q4 + pi/2;
theta5 = q1;
theta6 = 0;
theta7 = -pi/2;
theta8 = q4 + pi/2;

d1 = L4;
d2 = L5 + L6 + q2;
d3 = L7 + L8 + q3;
d4 = L10 + L11 + L13;
d5 = L4;
d6 = L5 + q2 + L6/2;
d7 = L7 + q3 + L8 + L14;
d8 = L10 + L11 +L13;

tz1 = sym('tz1',[4,4]);
tz1(3,4) = d1;
tz1 = tz1 + eye(4);

tz2 = sym('tz2',[4,4]);
tz2(3,4) = d2;
tz2 = tz2 + eye(4);


tz3 = sym('tz3',[4,4]);
tz3(3,4) = d3;
tz3 = tz3 + eye(4);

tz4 = sym('tz4',[4,4]);
tz4(3,4) = d4;
tz4 = tz4 + eye(4);

tz5 = sym('tz5',[4,4]);
tz5(3,4) = d5;
tz5 = tz5 + eye(4);

tz6 = sym('tz6',[4,4]);
tz6(3,4) = d6;
tz6 = tz6 + eye(4);

tz7 = sym('tz7',[4,4]);
tz7(3,4) = d7;
tz7 = tz7 + eye(4);

tz8 = sym('tz8',[4,4]);
tz8(3,4) = d8;
tz8 = tz8 + eye(4);



tx1 = sym('tx1',[4,4]);
tx1(1,4) = a1;
tx1 = tx1 + eye(4);

tx2 = sym('tx2',[4,4]);
tx2(1,4) = a2;
tx2 = tx2 + eye(4);


tx3 = sym('tx3',[4,4]);
tx3(1,4) = a3;
tx3 = tx3 + eye(4);

tx4 = sym('tx4',[4,4]);
tx4(1,4) = a4;
tx4 = tx4 + eye(4);

tx5 = sym('tx5',[4,4]);
tx5(1,4) = a5;
tx5 = tx5 + eye(4);

tx6 = sym('tx6',[4,4]);
tx6(1,4) = a6;
tx6 = tx6 + eye(4);

tx7 = sym('tx7',[4,4]);
tx7(1,4) = a7;
tx7 = tx7 + eye(4);

tx8 = sym('tx8',[4,4]);
tx8(1,4) = a8;
tx8 = tx8 + eye(4);



%% Joints and End-effector
% TODO: DEFINE THE Robot Base SEE FIG.6 
T0_W= [0 -1 0; 0 0 -1; 1 0 0 ];

% TODO: DEFINE THE Relative Homogeneous Transformations (symbolic form)

T1_0 = RotZ(theta1) * tz1 * tx1 + RotX(alpha1);
				
T2_1 = RotZ(theta2) * tz2 * tx2 + RotX(alpha2);
				
T3_2 = RotZ(theta3) * tz3 * tx3 + RotX(alpha3);

T4_3 = RotZ(theta4) * tz14 * tx4 + RotX(alpha4);






% TODO: DEFINE Homogeneous Transformations wrt BASE frame (Numeric computation)

T1_0 = T1_0;
T2_0 = (T1_0 * T2_1)^(-1) ;
T3_0 = (T1_0 * T2_1 * T3_2)^(-1);
T4_0 = (T1_0 * T2_1 * T3_2 * T4_3)^(-1);

% TODO: DEFINE Homogeneous Transformations wrt WORLD frame (Numeric computation)

T1_W = T0_W * T1_0;
T2_W = T0_W * T1_0 * T2_1;
T3_W = T0_W * T1_0 * T2_1 * T3_2;
T4_W = T0_W * T1_0 * T2_1 * T3_2 * T4_3;

%Compute the POSE of end-effector with respect to the world coordinate
%frame
% TODO
Xef_W=FK_robot([q1;q2;q3;q4;L1;L2;Ln]);

%DRAW Base, Links, CFs and EF
%Use the functions from tutorial 1 to plot, joints, links, all CFs (base, end-effector) 


% %% Centers of Mass
% % TODO: Relative Homogeneous Transformations for each CM (symbolic equations)
% Tcm1_0=??;		
% Tcm2_1=??;				
% Tcm3_2=??;
% Tcm4_3=??;
% 
% % TODO: Homogeneous Transformations wrt base frame (Numeric computation)
% Tcm1_0=eye(4);
% Tcm2_0=eye(4);
% Tcm3_0=eye(4);
% Tcm4_0=eye(4);
% 
% 
% % TODO: Homogeneous Transformations wrt World frame (Numeric computation)
% Tcm1_W=eye(4);
% Tcm2_W=eye(4);
% Tcm3_W=eye(4);
% Tcm4_W=eye(4);
% 
% 
% %DRAW CFs CMs
% 
% %TODO: Use the functions from tutoria 1 to plot CF for each center of mass
% %(3 in total). Create a matlab function called plotcfs.m
% 
% 
% %TODO: Compute the POSE of the cm 2 with respect to the robot base (CF 0) [numeric value]
% Xcm2_0=zeros(6,1);
% 
% % TODO Output ONLY the position vector for end effector and cm_2 (output size 6X1) 
% Xout=[Xef;Xcm2];




end

