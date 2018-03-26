clear
clc

%% Symbols
syms q1 q2 q3 qp1 qp2 qp3 qpp1 qpp2 qpp3 L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 m1 m2 m3 g gx gy gz real
syms qr1 qr2 qr3 qrp1 qrp2 qrp3 qrpp1 qrpp2 qrpp3 real
syms I111 I112 I113 I122 I123 I133 real
syms I211 I212 I213 I222 I223 I233 real
syms I311 I312 I313 I322 I323 I333 real

% L11 = sym(0.015);
L11 = sym(0.1157);
%% Denevit-Hartenberg-Table
%% Old DH Table
% DH = [...
%            q1,        L1,        0, -pi/2;...
%       q2-pi/2,         0,       L3,    pi;...
%           -q3,  -(L2+L4),       L5,    pi;...
%             0,        L6,        0,     0;...
%       q2-pi/2,        L7,       L8,     0;...
%            q3, -(L9+L11),      L10,     0;...
%      ];

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

H1_0 = HT{2,1};
H2_0 = HT{2,2};
H3_0 = HT{2,3};

Hcm1_0 = HT{4,1};
Hcm2_0 = HT{4,2};
Hcm3_0 = HT{4,3};


%% Plot Robot

A = cell2sym(HT);
T = subs(A,[q1 q2 q3 L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11],[0,0,0,0.128,0.1639,0.6127,0.0922,0.6873,0.1,0.15,0.6127/2,0.01,0.3437,0.1157]);
A = mat2cell(T,[4,4,4,4],[4,4,4]);

figure(1);
clf
view(5,15);
robotPlot(A);
%% Differential Kinematics

tcm1_0 = HT{4,1}(1:3,4);

tcm2_0 = HT{4,2}(1:3,4);
tcm2_1 = tcm2_0 - HT{2,1}(1:3,4);

tcm3_0 = HT{4,3}(1:3,4);
tcm3_1 = tcm3_0 - HT{2,1}(1:3,4);
tcm3_2 = tcm3_0 - HT{2,2}(1:3,4);

tef_0 = HT{2,3}(1:3,4);
tef_1 = tef_0 - HT{2,1}(1:3,4);
tef_2 = tef_0 - HT{2,2}(1:3,4);

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
                
Jef = [...
cross( z0,tef_0),  cross(z1,tef_1),  cross(z2,tef_2);...
              z0,               z1,              z2];
          
Jefd = diff(Jef,q1)*qp1 + diff(Jef,q2)*qp2 + diff(Jef,q3)*qp3;

%% Dynamic Model (M(q)qdd + C(q,qd)qd + G(q) = tau)

                 %%%%%%%%%%%%%%% M(q)%%%%%%%%%%%%%%%%
M = sym(zeros(3));

m = [m1,m2,m3];
J = {Jcm1,Jcm2,Jcm3};
I1 = [I111, I112, I113; I112, I122, I123; I113, I123, I133];
I2 = [I211, I212, I213; I212, I222, I223; I213, I223, I233];
I3 = [I311, I312, I313; I312, I322, I323; I313, I323, I333];
I = {I1, I2, I3};

Rcm1 = HT{4,1}(1:3,1:3);
Rcm2 = HT{4,2}(1:3,1:3);
Rcm3 = HT{4,3}(1:3,1:3);
R = {Rcm1, Rcm2, Rcm3};

for i = 1:3
    M = M + (m(i)*J{i}(1:3,:)'*J{i}(1:3,:)) + (J{i}(4:6,:)'*R{i}*I{i}*R{i}'*J{i}(4:6,:));
end 
M = simplify(M);



               %%%%%%%%%%%%%%% C(q,qp)%%%%%%%%%%%%%%%%
C = sym('C', [3,3],'real');
q = [q1, q2, q3];
qd = [qp1, qp2, qp3];
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


Md = diff(M,q1)*qp1 + diff(M,q2)*qp2 + diff(M,q3)*qp3;
N = simplify(Md - 2*C);
x = sym('x', [3,1],'real');
test = simplify(x'*N*x);


               %%%%%%%%%%%%%%%%% G(q)%%%%%%%%%%%%%%%%%%%%
               
g_vector = g*[gx;gy;gz];
tcm = [tcm1_0, tcm2_0, tcm3_0];

P = sym(0);
for i = 1:3
    P = P + (m(i)*g_vector'*tcm(:,i));
end

G = [diff(P,q1); diff(P,q2); diff(P,q3)];




%% Dynamic equations

qpp = [qpp1;qpp2;qpp3];
qp = [qp1; qp2; qp3];

qppr = [qrpp1; qrpp2; qrpp3];
qpr = [qrp1; qrp2; qrp3];

exp = M*qpp + C*qp + G;
expr = M*qppr + C*qpr + G;


%% Regressor

params = [L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,m1,m2,m3,g,gx,gy,gz,I111,I112,I113,I122,I123,I133,I211,I212,I213,I222,I223,I233,I311,I312,I313,I322,I323,I333];


[Y,Theta] = regressor(exp, params);
[Yr,Thetar] = regressor(expr, params);

test1 = simplify(exp - Y*Theta);
test2 = simplify(expr - Yr*Thetar);


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
    





