
clear
clc

%% Symbols
syms q1 q2 q3 q1p q2p q3p q1pp q2pp q3pp L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11 m1 m2 m3 g real
syms I11_1 I11_2 I11_3 I12_2 I12_3 I13_3 real
syms I21_1 I21_2 I21_3 I22_2 I22_3 I23_3 real
syms I31_1 I31_2 I31_3 I32_2 I32_3 I33_3 real


%% Denevit-Hartenberg-Table
DH = [...
           q1,        L1,        0,  pi/2;...
        q2+pi,         0,       L3,     0;...
           q3,     L2+L4,       L5,     0;...
           q1,        L6,        0,     0;...
        q2+pi,        L7,       L8,     0;...
           q3,        L9,      L10,     0;...
     ];
 


 
%% Forward Kinematics
HT = DH2H(DH);

%%
H1_0 = subs(HT{2,1},[L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11],[0.1273;0.1639;0.6127;0.0922;0.6873;0.1;0.15;0.6127/2;0.01;0.3437;0.15]');
H2_0 = subs(HT{2,2},[L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11],[0.1273;0.1639;0.6127;0.0922;0.6873;0.1;0.15;0.6127/2;0.01;0.3437;0.15]');
H3_0 = subs(HT{2,3},[L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11],[0.1273;0.1639;0.6127;0.0922;0.6873;0.1;0.15;0.6127/2;0.01;0.3437;0.15]');

Hcm1_0 = subs(HT{4,1},[L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11],[0.1273;0.1639;0.6127;0.0922;0.6873;0.1;0.15;0.6127/2;0.01;0.3437;0.15]');
Hcm2_0 = subs(HT{4,2},[L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11],[0.1273;0.1639;0.6127;0.0922;0.6873;0.1;0.15;0.6127/2;0.01;0.3437;0.15]');
Hcm3_0 = subs(HT{4,3},[L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11],[0.1273;0.1639;0.6127;0.0922;0.6873;0.1;0.15;0.6127/2;0.01;0.3437;0.15]');;


q = [q1 q2 q3];
matlabFunction(H1_0, H2_0, H3_0, Hcm1_0, Hcm2_0, Hcm3_0,...
                   'File','../H', ...
                   'Vars',{q}, ...
                   'Optimize',false, ...
                   'Comments','Compute Transformation Matrices: H(q)', ...
                   'Outputs',{'H1_0','H2_0','H3_0', 'Hcm1_0','Hcm2_0','Hcm3_0'});

%% Differential Kinematics
 
tcm1_0 = HT{4,1}(1:3,4);

tcm2_0 = HT{4,2}(1:3,4);
tcm2_1 = tcm2_0 - HT{2,1}(1:3,4);

tcm3_0 = HT{4,3}(1:3,4);
tcm3_1 = tcm3_0 - HT{2,1}(1:3,4);
tcm3_2 = tcm3_0 - HT{2,2}(1:3,4);

z0 = [0;0;1];
z1 = HT{2,1}(1:3,3);
z2 = HT{2,2}(1:3,3);

Jcm1 = [...
      cross(z0,tcm1_0), zeros(3,1),  zeros(3,1);...
                    z0, zeros(3,1),  zeros(3,1)];     
               
Jcm2 = [...
      cross(z0,tcm2_0), cross(z1,tcm2_1),  zeros(3,1);...
                    z0,               z1,  zeros(3,1)]; 
              
Jcm3 = [...
      cross(z0,tcm3_0), cross(z1,tcm3_1), cross(z2,tcm3_2);...
                    z0,               z1,              z2];



%% Dynamic Model (M(q)qdd + C(q,qd)qd + G(q) = tau)

                 %%%%%%%%%%%%%%% M(q)%%%%%%%%%%%%%%%%
M = sym(zeros(3));

m = [m1,m2,m3];
J = {Jcm1,Jcm2,Jcm3};
I1 = [I11_1, I11_2, I11_3; I11_2, I12_2, I12_3; I11_3, I12_3, I13_3];
I2 = [I21_1, I21_2, I21_3; I21_2, I22_2, I22_3; I21_3, I22_3, I23_3];
I3 = [I31_1, I31_2, I31_3; I31_2, I32_2, I32_3; I31_3, I32_3, I33_3];
I = {I1, I2, I3};

Rcm1 = HT{4,1}(1:3,1:3);
Rcm2 = HT{4,2}(1:3,1:3);
Rcm3 = HT{4,3}(1:3,1:3);
R = {Rcm1, Rcm2, Rcm3};

for i = 1:3
    M = M + (m(i)*J{i}(1:3,:)'*J{i}(1:3,:)) + (J{i}(4:6,:)'*R{i}*I{i}*R{i}'*J{i}(4:6,:));
end 
M = simplify(M);



               %%%%%%%%%%%%%%% C(q,qd)%%%%%%%%%%%%%%%%
C = sym('C', [3,3],'real');
q = [q1, q2, q3];
qd = [q1p, q2p, q3p];
temp = sym('temp', [1,3],'real');
for k = 1:3
    for j = 1:3
        for i = 1:3
            temp(i) = (diff(M(k,j),q(i)) + diff(M(k,i),q(j)) - diff(M(i,j),q(k)))*qd(i);
        end
        C(k,j) = 1/2 * sum(temp);
    end
end
C = simplify(C);


Md = diff(M,q1)*q1p + diff(M,q2)*q2p + diff(M,q3)*q3p;
N = simplify(Md - 2*C);
x = sym('x', [3,1],'real');
test = simplify(x'*N*x);


               %%%%%%%%%%%%%%%%% G(q)%%%%%%%%%%%%%%%%%%%%
               
g_vector = [0;0;g];
tcm = [tcm1_0, tcm2_0, tcm3_0];

P = sym(0);
for i = 1:3
    P = P + (m(i)*g_vector'*tcm(:,i));
end

G = [diff(P,q1); diff(P,q2); diff(P,q3)];

%%

q = [q1; q2; q3]; 
qp = [q1p; q2p; q3p];

L = [L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11];
m = [m1 m2 m3];
I1 = [I11_1 I11_2 I11_3 I12_2 I12_3 I13_3];
I2 = [I21_1 I21_2 I21_3 I22_2 I22_3 I23_3];
I3 = [I31_1 I31_2 I31_3 I32_2 I32_3 I33_3];

L_num = [0.1273;0.1639;0.6127;0.0922;0.6873;0.1;0.15;0.6127/2;0.01;0.3437;0.15]';
m_num = [5 15 15];
I1_num = [0.04;0.02;0.02;0.04;0.02;0.04]';
I2_num = I1_num;
I3_num = I2_num;
g_num = 9.81;


M_num = subs(M, [L,m,I1,I2,I3,g], [L_num,m_num,I1_num,I2_num,I3_num,g_num]);
C_num = subs(C, [L,m,I1,I2,I3,g], [L_num,m_num,I1_num,I2_num,I3_num,g_num]);
G_num = subs(G, [L,m,I1,I2,I3,g], [L_num,m_num,I1_num,I2_num,I3_num,g_num]);

matlabFunction(M_num, C_num, G_num,...
                   'File','../EulerLagrangeMatrices', ...
                   'Vars',{q,qp}, ...
                   'Outputs',{'M','C','G'});


%% Get Final Matrix Equation

qpp = [q1pp;q2pp;q3pp];
qp = [q1p; q2p; q3p];

tau = M*qpp + C*qp +G;


%% Regressor
params = [L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,m1,m2,m3,g,I11_1,I11_2,I11_3,I12_2,I12_3,I13_3,I21_1,I21_2,I21_3,I22_2,I22_3,I23_3,I31_1,I31_2,I31_3,I32_2,I32_3,I33_3];



[Y,Theta] = regressor(tau, params);

% Validate regressor
test1 = simplify(tau - Y*Theta);

if(norm(test1) == 0)
    fprintf('Found Regressor!\n');
else
    fprintf('Error in Regressor!\n')
end


clearvars

%% Functions
function [Y,Theta] = regressor(exp,params)
    [c1,t1] = coeffs(exp(1),params);
    [c2,t2] = coeffs(exp(2),params);
    [c3,t3] = coeffs(exp(3),params);

    [t2,c2] = combine(t1,t2,c2);
    idx=find(t2 == 0);
    t2(idx) = t1(idx);
    [t3,c3] = combine(t2,t3,c3);
    idx=find(t3 == 0);
    t3(idx) = t2(idx);
    [~,c1] = combine(t3,t1,c1);

    Theta = t3';
    Y = [c1;c2;c3];
end

function [t2,c2] = combine(t1,t2,c2)
    for i = 1:size(t1,2)
        idx = find(t2 == t1(i));
         % reorder t2
        if(isempty(idx))
            % insert zero
            t2 = [t2(1:i-1),0,t2(i:end)];
            c2 = [c2(1:i-1),0,c2(i:end)];
        elseif(idx ~= i)
            % swap if in wrong position
            t2([idx i]) = t2([i idx]);
            c2([idx i]) = c2([i idx]);
        end
    end
end


