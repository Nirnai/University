clear
clc

%% Symbols
syms q1 q2 q3 qp1 qp2 qp3 qpp1 qpp2 qpp3 L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 m1 m2 m3 g gx gy gz real
syms I111 I112 I113 I122 I123 I133 real
syms I211 I212 I213 I222 I223 I233 real
syms I311 I312 I313 I322 I323 I333 real

L11 = sym(0.015);
%% Denevit-Hartenberg-Table
DH = [...
           q1,        L1,        0, -pi/2;...
      q2-pi/2,         0,       L3,    pi;...
          -q3,  -(L2+L4),       L5,    pi;...
            0,        L6,        0,     0;...
       q2-pi/2,        L7,       L8,     0;...
           q3, -(L9+L11),      L10,     0;...
     ];


 
%% Forward Kinematics
HT = DH2H(DH);

H1_0 = HT{2,1};
H2_0 = HT{2,2};
H3_0 = HT{2,3};

Hcm1_0 = HT{4,1};
Hcm2_0 = HT{4,2};
Hcm3_0 = HT{4,3};


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

%% Plot Robot

A = cell2sym(HT);
T = subs(A,[q1 q2 q3 L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11],[0,0,0,0.128,0.1639,0.6127,0.0922,0.6873,0.1,0.15,0.6127/2,0.01,0.3437,0.15]);
A = mat2cell(T,[4,4,4,4],[4,4,4]);

figure(1);
clf
view(5,15);
robotPlot(A);

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




% %% Get Final Matrix Equation
% 
% qdd = [qdd1;qdd2;qdd3];
% qd = [qd1; qd2; qd3];
% 
% exp = M*qdd + C*qd +G;
% 
% %% Create a Parameter Vector
% 
% params = [L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11,m1,m2,m3,g,I11_1,I11_2,I11_3,I12_2,I12_3,I13_3,I21_1,I21_2,I21_3,I22_2,I22_3,I23_3,I31_1,I31_2,I31_3,I32_2,I32_3,I33_3];



% %% find Regressor
% param = cell(1,5);
% regressor = cell(1,3);
% 
% for i = 1:3
%     [regressor{1,i},param{1,i}] = coeffs(exp(i),params);    
% end
% 
% % %%%%%%%%%
% % param_test = cell(1,3);
% % regressor_test = cell(1,3);
% % for i = 1:3
% %     [regressor_test{1,i},param_test{1,i}] = coeffs(exp(i),params);    
% % end
% % %%%%%%%%%%%
% 
% param{1,4} = param{1,2};
% 
% for i = 1:size(param{1,1},2)
%     idx = find(param{1,2} == param{1,1}(1,i));
%     if(isempty(idx))
%         param{1,2} = [param{1,2}(1:i-1),0,param{1,2}(i:end)];
%         regressor{1,2} = [regressor{1,2}(1:i-1),0,regressor{1,2}(i:end)];
%         param{1,4} = [param{1,4}(1:i-1),param{1,1}(1,i),param{1,4}(i:end)];        
%     else
%         param{1,2}([i idx]) = param{1,2}([idx i]);
%         regressor{1,2}([i idx]) = regressor{1,2}([idx i]);
%         param{1,4}([i idx]) = param{1,4}([idx i]);
%     end
% end
% 
% 
% param{1,5} = param{1,3};
% 
% for i = 1:size(param{1,4},2)
%     idx = find(param{1,3} == param{1,4}(1,i));
%     if(isempty(idx))
%         param{1,3} = [param{1,3}(1:i-1),0,param{1,3}(i:end)]; 
%         regressor{1,3} = [regressor{1,3}(1:i-1),0,regressor{1,3}(i:end)];
%         param{1,5} = [param{1,5}(1:i-1),param{1,4}(1,i),param{1,5}(i:end)];
%     elseif(i == idx)
%     else
%         param{1,3}([i idx]) = param{1,3}([idx i]);
%         regressor{1,3}([i idx]) = regressor{1,3}([idx i]);
%         param{1,5}([i idx]) = param{1,5}([idx i]);
%     end
% end
%     
% 
% for i = 1:size(param{1,5},2)
%     idx = find(param{1,1} == param{1,5}(1,i));
%     if(isempty(idx))
%         param{1,1} = [param{1,1}(1:i-1),0,param{1,1}(i:end)]; 
%         regressor{1,1} = [regressor{1,1}(1:i-1),0,regressor{1,1}(i:end)];
%     else
%         param{1,1}([i idx]) = param{1,1}([idx i]);
%         regressor{1,1}([i idx]) = regressor{1,1}([idx i]);
%     end
% end
% 
% 
% Y = sym('Y', [3,size(regressor{1,1},2)]);
% 
% Y(1,:) = regressor{1,1};
% Y(2,:) = regressor{1,2};
% Y(3,:) = regressor{1,3};
% 
% Theta = param{1,5}';
% 
