%% RPPR Symbols
syms qd1 qd2 qd3 qd4 q1 q2 q3 q4 
% syms L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11 L12 L13 L14 L15 L16 L17 L18
% syms c T1 T2 T3

global L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11 L12 L13 L14 L15 L16 L17 L18
global c T1 T2 T3


% Load Robot Parameters into global variables
parameters(0);

%% D-H-Table Joints

DH = [...
                q1,            0,      0,      -pi/2;...
             -pi/2,     L5+L6+q2,      0,       pi/2;...
                 0, L13+L7+L8+q3,    L14,  (T1+pi/2);...
          -pi/2+q4,  L10+L11+L15,   -L16,          0; ... 
           pi/2+q1,            0,  -L5/2,          0;...
                 0, L5+(L6/2)+q2,   L7/2,          0;...
           pi/2+T3, L7+L8+L13+q3,    L18,          0;...
        -(pi/2-q4),  L10+L11+L15,   -L16,          0 ...
        ];
    
       
%% Forward Kinematics
H0_W= [ 0 -1  0      L1;...
        0  0 -1 (L2-L4); ...
        1  0  0      L3; ...
        0  0  0      1];
% Compute Symbilic Transformations
HT = DH2H(DH);

H1_W   = H0_W * HT{2,1};
H2_W   = H0_W * HT{2,2};
H3_W   = H0_W * HT{2,3};
H4_W   = H0_W * HT{2,4};


T = sym(zeros(4,4,5));
T(:,:,1) = H0_W;
T(:,:,2) = H1_W;
T(:,:,3) = H2_W;
T(:,:,4) = H3_W;
T(:,:,5) = H4_W;

T = subs(T,[q1,q2,q3,q4],[0,0,0,0]);

view(5,15);
robotPlot(T);



% Generate Function to Evaluate Transformations during Runtime 
q = [q1; q2; q3; q4];
matlabFunction(H0_W, H1_W, H2_W, H3_W, H4_W, ...
                   'File','../H', ...
                   'Vars',{q}, ...
                   'Optimize',false, ...
                   'Comments','Compute Transformation Matrices: H(q)', ...
                   'Outputs',{'H0_W','H1_W','H2_W','H3_W','H4_W'});


%% Forward Kinematics
% Endeffector Pose
t_ef_0 = simplify(HT{2,4}(1:3,4));
[r_ef_0,p_ef_0,y_ef_0] = rot2eul(HT{2,4}(1:3,1:3));
Xef = [t_ef_0;r_ef_0;p_ef_0;y_ef_0];

% CenterOfMass 2 Pose
t_cm2_0 = simplify(HT{4,2}(1:3,4));
[r_cm2_0,p_cm2_0,y_cm2_0] = rot2eul(HT{4,2}(1:3,1:3));
Xcm2 = [t_cm2_0;r_cm2_0;p_cm2_0;y_cm2_0];




% Create Forward Kinematics Function
matlabFunction(Xef, 'File', '../f', 'Vars', {q},'Comments',' ForwarKinematics x = f(q)', 'Optimize',false);
matlabFunction(Xcm2, 'File', '../fcm2', 'Vars', {q},'Comments',' ForwarKinematics xcm2 = fcm2(q)', 'Optimize',false);


%% Differential Kinematics
t_ef_0 = HT{2,4}(1:3,4);
t_ef_1 = t_ef_0 - HT{2,1}(1:3,4);
t_ef_2 = t_ef_0 - HT{2,2}(1:3,4);
t_ef_3 = t_ef_0 - HT{2,3}(1:3,4);

t_cm1_0 = HT{4,1}(1:3,4);

t_cm2_0 = HT{4,2}(1:3,4);
t_cm2_1 = HT{4,2}(1:3,4);

t_cm3_0 = HT{4,3}(1:3,4);
t_cm3_1 = t_cm3_0 - HT{2,1}(1:3,4);
t_cm3_2 = t_cm3_0 - HT{2,2}(1:3,4);

z0 = [0;0;1];
z1 = HT{2,1}(1:3,3);
z2 = HT{2,2}(1:3,3);
z3 = HT{2,3}(1:3,3);

% Jacobi Matrix 
% Base - Endeffector
J_ef_0 = [...
            cross(z0,t_ef_0),         z1,         z2, cross(z3,t_ef_3);
                          z0, zeros(3,1), zeros(3,1),              z3];
              
% Base - CenterOfMass1
J_cm1_0 = [...
    cross(z0,t_cm1_0), zeros(3,1), zeros(3,1), zeros(3,1);
                   z0, zeros(3,1), zeros(3,1), zeros(3,1)];              
              
% Base - CenterOfMass2
J_cm2_0 = [...
    cross(z0,t_cm2_0),          z1, zeros(3,1), zeros(3,1);
                   z0,  zeros(3,1), zeros(3,1), zeros(3,1)];              
             
% Base - CenterOfMass3
J_cm3_0 = [...
    cross(z0,t_cm3_0),          z1,          z2, zeros(3,1);
                   z0,  zeros(3,1),  zeros(3,1), zeros(3,1)];


matlabFunction(J_ef_0, 'File', '../J', 'Vars', {q},'Comments',' Jacobian: J(q)', 'Optimize',false);

matlabFunction(J_cm2_0, 'File', '../Jcm2', 'Vars', {q}, 'Comments',' Jacobian: Jcm2(q)', 'Optimize',false);

clearvars
               