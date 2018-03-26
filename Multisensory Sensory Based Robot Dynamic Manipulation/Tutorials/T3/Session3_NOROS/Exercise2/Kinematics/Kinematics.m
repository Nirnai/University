%% Variables and Parameters
syms q1p q2p q3p q1 q2 q3 
global L1 L2 L3 L4 L5 L6

parameters(0);

q  = [q1;q2;q3];
qp = [q1p;q2p;q3p];

%% Denevit-Hartenberg Table

DH =[...
           q1,         0,        0,  pi/2;...
      q2+pi/2,         0,       L2,    pi;...
           q3,  -(L3+L5),    L4+L6,    pi;...
    ];


%% Homogeneous Transdormations

HT = DH2H(DH);

% Define the Robot Base
H0_W=[ 1, 0, 0,  0; ...
       0, 1, 0,  0; ...
       0, 0, 1, L1; ...
       0, 0, 0,  1];
   
H1_W   = H0_W * HT{2,1};
H2_W   = H0_W * HT{2,2};
H3_W   = H0_W * HT{2,3};

matlabFunction(H0_W, H1_W, H2_W, H3_W, ...
                   'File','../H_UR10_3', ...
                   'Vars',{q}, ...
                   'Optimize',false, ...
                   'Comments','Compute Transformation Matrices: H(q)', ...
                   'Outputs',{'H0_W','H1_W','H2_W','H3_W'});



%% Forward Kinematics - End Effector
% Forward Kinematics
t = HT{2,3}(1:3,4);
[r,p,y] = rot2eul(HT{2,3}(1:3,1:3));
Xef = simplify([t;r;p;y]);

% Jacobian
t_ef_0 = HT{2,3}(1:3,4);
t_ef_1 = t_ef_0 - HT{2,1}(1:3,4);
t_ef_2 = t_ef_0 - HT{2,2}(1:3,4);

z0 = [0;0;1];
z1 = HT{2,1}(1:3,3);
z2 = HT{2,2}(1:3,3);    

Jef = simplify([ ...
        cross(z0,t_ef_0), cross(z1, t_ef_1), cross(z2, t_ef_2);
                      z0,                z1,                z2
      ]);


  
Xefp =  simplify(Jef * qp);

matlabFunction(Xef, 'File','../FK_UR10_3',   'Vars', {q},'Comments',' x = f(q)', 'Optimize',false);
matlabFunction(Jef, 'File','../Jef_0_UR10_3','Vars', {q},'Comments',' J(q)', 'Optimize',false);
matlabFunction(Xefp,'File','../DifFK_UR10_3','Vars', {q,qp},'Comments',' Xp = J(q) * qp', 'Optimize',false);
  
clearvars