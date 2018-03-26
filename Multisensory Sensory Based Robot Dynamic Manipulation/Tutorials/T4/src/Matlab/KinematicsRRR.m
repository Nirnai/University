addpath('Subfunctions/');

%%

syms q1 q2 q3 qd1 qd2 qd3 q1dd q2dd q3dd m1 m2 m3 l1 l2 l3 k1 k2 k3 z0 T1 T2 b1 b2 b3 g real

%% D-H-Table
DH = [...
          0, l1+q1,    0, -pi/2;...
    q2+pi/2,     0,  -l2,     0;...
         q3,     0,  -l3,     0;...
          0,    q1,    0,     0;...
    q2+pi/2,     0,  -l2,     0;...
         q3,     0,   l3,     0;...
    ];


%%

% Compute Symbilic Transformations
HT = DH2H(DH);

H0_W = [...
        1  0  0 0;
        0 -1  0 0;
        0  0 -1 0;
        0  0  0 1];
   
H1_W = H0_W * HT{2,1};
H2_W = H0_W * HT{2,2};
H3_W = H0_W * HT{2,3};


% Generate Function to Evaluate Transformations during Runtime 
u = [q1 q2 q3 qd1 qd2 qd3 l1 l2 l3 m1 m2 m3 k1 k2 k3 z0 T1 T2 b1 b2 b3 g];
matlabFunction(H1_W, H2_W, H3_W,...
                   'File','HTs', ...
                   'Vars',{u}, ...
                   'Outputs',{'H1_W','H2_W','H3_W'});
               

               
%%
%% Plot Robot
% 
% A = cell2sym(HT);
% T = subs(A,[q1 q2 q3 l1 l2 l3 ],[0,0,0,1.575,1.075,1.075]);
% A = mat2cell(T,[4,4,4,4],[4,4,4]);
% 
% figure(1);
% clf
% view(5,15);
% robotPlot(A);

% %%
