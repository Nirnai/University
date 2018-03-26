addpath('Subfunctions/');

%% RPPR Symbols
syms qd1 qd2 qd3 qd4 q1 q2 q3 q4 L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11 L12 c T1 T2

L13 = L9 * sin(T1);
L14 = L9 * cos(T1);
L15 = (c/2 + L12) * cos(T1);
L16 = (c/2 + L12) * sin(T1);
L17 = L10 * 1/cos(T1);
T3 = atan(L14/L17);
L18 = L17 * 1/cos(T3);

%% D-H-Table Joints
DH_j = [...
                q1,            0,    0,      -pi/2;...
             -pi/2,     L5+L6+q2,    0,       pi/2;...
                 0, L13+L7+L8+q3,  L14,  (T1+pi/2);...
        -pi/2+q4,  L10+L11+L15,   -L16,          0 ... 
       ];

%% D-H-Table Center of Masses
DH_cm = [...
            pi/2+q1,            0, -L5/2,   0;...
                  0, L5+(L6/2)+q2,  L7/2,   0;...
            pi/2+T3, L7+L8+L13+q3,   L18,   0;...
         -(pi/2-q4),  L10+L11+L15,   -L16,   0 ...
        ];
    
%% Forward Kinematics
H0_W= [ 0 -1  0      L1;...
        0  0 -1 (L2-L4); ...
        1  0  0      L3; ...
        0  0  0      1];
% Compute Symbilic Transformations
[HT_Joints, HT_CM] = DH2H(DH_j,DH_cm,[q1 q2 q3 q4 L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11 L12 c T1 T2]);

% Generate Function to Evaluate Transformations during Runtime 
u = [qd1; qd2; qd3; qd4; q1; q2; q3; q4; L1; L2; L3; L4; L5; L6; L7; L8; L9; L10; L11; L12; c; T1; T2];
matlabFunction(H0_W, HT_Joints{2,1}, HT_Joints{2,2}, HT_Joints{2,3}, HT_Joints{2,4}, HT_CM{2,2}, ...
                   'File','Subfunctions/HTs', ...
                   'Vars',{u}, ...
                   'Outputs',{'H0_W','H1_0','H2_0','H3_0','H4_0', 'Hcm2_0'});


%% Endeffector Pose

t = simplify(HT_Joints{2,4}(1:3,4));
[r,p,y] = rot2eul(HT_Joints{2,4}(1:3,1:3));

x_ef = [t;r;p;y];

