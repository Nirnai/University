clearvars

syms L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11 L12 real;
syms q1 q2 q3 q4 q5 q6 real;
syms q1p q2p q3p q4p q5p q6p t real;

% L = [L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11 L12];
% q = [q1,q2,q3,q4,q5,q6];
% qip = [q1p,q2p,q3p,q4p,q5p,q6p];
% 
% numDOF = 6;
% 
% % DH(1,:) = [q1, L1, 0, pi/2];
% % DH(2,:) = [q2+pi, 0, L3, 0];
% % DH(3,:) = [q3+pi/2, L2, L5, 0];
% % DH(4,:) = [q4, 0, 0, -pi/2];
% % DH(5,:) = [q5, -L11, 0, pi/2];
% % DH(6,:) = [q6, L12, 0, 0];
% 
% DH(1,:) = [q1, L1, 0, pi/2];
% DH(2,:) = [q2+pi/2, 0, L3, 0];
% DH(3,:) = [q3, 0, L5, 0];
% DH(4,:) = [q4-pi/2, L2, 0, -pi/2];
% DH(5,:) = [q5, L11, 0, pi/2];
% DH(6,:) = [q6, L12+L4, 0, 0];
% 
% for i = 1:numDOF
%    
%     if has(DH(i,1),q(i))
%         tagDOF(i,1) = 'r';
%     else
%         tagDOF(i,1) = 'p';
%     end
%         
% end
% 
% % DH_cm(1,:) = [q1, L6, 0, 0];
% % DH_cm(2,:) = [q2+pi, L7, L8, 0];
% % DH_cm(3,:) = [q3, L9, L10, 0];
% % DH_cm(4,:) = [q4, 0, 0, 0];
% % DH_cm(5,:) = [q5, -L11, 0, 0];
% % DH_cm(6,:) = [q6, L12, 0, 0];
% 
% DH_cm(1,:) = [q1, L6, 0, 0];
% DH_cm(2,:) = [q2+pi/2, L7, L8, 0];
% DH_cm(3,:) = [q3, L9, L10, 0];
% DH_cm(4,:) = [q4, L2/2, 0, 0];
% DH_cm(5,:) = [q5, L11/2, 0, 0];
% DH_cm(6,:) = [q6, (L4+L12)/2, 0, 0];
% 
% % T0_W=sym([1,0,0,-0.5;
% % 0,1,0,-0.5;
% % 0,0,1,0;
% % 0,0,0,1]);
% T0_W=sym(eye(4));
% 
% T0_W = simplify(T0_W);
% % DEFINE THE Relative Homogeneous Transformations 
% for i = 1:numDOF 
%     Hz = [RotZ(DH(i,1)),[0,0,0]';[0,0,0,1]];
%     Tz(:,4) = [0, 0, DH(i,2), 0]';
%     Tz = Tz +eye(4);
%     Tx(:,4) = [DH(i,3), 0, 0, 0]';
%     Tx = Tx + eye(4);
%     Hx = [RotX(DH(i,4)),[0,0,0]';[0,0,0,1]];
%     
%     Hi_iminus1(:,:,i) = Hz*Tz*Tx*Hx;
%     Tz = Tz -eye(4);
%     Tx = Tx - eye(4);
% end
% Hi_iminus1 = simplify(Hi_iminus1);
% 
% 
% % TODO: DEFINE Homogeneous Transformations wrt BASE frame (Numeric computation)
% 
% for i = 1:numDOF
%     
%    Hi_0(:,:,i) =  Hi_iminus1(:,:,i);
%    
%    j = i-1;
%    while j > 0
%        Hi_0(:,:,i) =  Hi_iminus1(:,:,j)*Hi_0(:,:,i);
%        j = j-1;
%    end
% end
% 
% Hi_0 = simplify(Hi_0);
% % TODO: DEFINE Homogeneous Transformations wrt WORLD frame (Numeric computation)
% 
% for i = 1:numDOF
%    
%     Hi_w(:,:,i) = T0_W*Hi_0(:,:,i);
% end
% Hi_w = simplify(Hi_w);
% 
% for i = 1:numDOF 
%     Hz = [RotZ(DH_cm(i,1)),[0,0,0]';[0,0,0,1]];
%     Tz(:,4) = [0, 0, DH_cm(i,2), 0];
%     Tz = Tz +eye(4);
%     Tx(:,4) = [DH_cm(i,3), 0, 0, 0];
%     Tx = Tx + eye(4);
%     Hx = [RotX(DH_cm(i,4)),[0,0,0]';[0,0,0,1]];
% 
%     Hcmi_iminus1(:,:,i) = Hz*Tz*Tx*Hx;
%     Tz = Tz -eye(4);
%     Tx = Tx - eye(4);
% end
% Hcmi_iminus1 = simplify(Hcmi_iminus1);
% 
% % TODO: Homogeneous Transformations wrt base frame (Numeric computation)
% for i = 1:numDOF
%     
%    Hcmi_0(:,:,i) =  Hcmi_iminus1(:,:,i);
%    
%    j = i-1;
%    while j > 0
%        Hcmi_0(:,:,i) =  Hi_iminus1(:,:,j)*Hcmi_0(:,:,i);
%        j = j-1;
%    end
% end
% Hcmi_0 = simplify(Hcmi_0);
% 
% 
% % TODO: Homogeneous Transformations wrt World frame (Numeric computation)
% 
% for i = 1:numDOF
%    
%     Hcmi_w(:,:,i) = T0_W*Hcmi_0(:,:,i);
% end
% 
% Hcmi_w = simplify(Hcmi_w);
% 
% 
% fid = fopen('MymatrixUR10Project.txt','wt');
% WriteToFile(fid,T0_W,'T0_W');
% WriteToFile(fid,Hi_iminus1,'Hi_iminus1');
% WriteToFile(fid,Hi_0,'Hi_0');
% WriteToFile(fid,Hi_w,'Hi_w');
% %Center of masses
% WriteToFile(fid,Hcmi_iminus1,'Hcmi_iminus1');
% WriteToFile(fid,Hcmi_0,'Hcmi_0');
% WriteToFile(fid,Hcmi_w,'Hcmi_w');
% fclose(fid);
% % TODO: Relative Homogeneous Transformations for each CM (symbolic equations)
% 
% 
% 
% %%  Define axis and radius wrt. world
% 
% kz = [0;0;1];
% ziminus1_w = sym(zeros(3,numDOF));
% ti_w = sym(zeros(3,numDOF));
% 
% ziminus1_w(:,1) = T0_W(1:3,1:3)*kz;
% for i = 2:numDOF
%    ziminus1_w(:,i) = Hi_w(1:3,1:3,i-1)*kz;
% end
% 
% for i = 1:numDOF
%    ti_w(:,i) = Hi_w(1:3,4,i);
% end
% 
% 
% %% Calculate End Effector Jacobian
% 
% Jv = sym(zeros(3,numDOF));
% Jw = sym(zeros(3,numDOF));
% 
% for i = 1:numDOF
%     if strcmp(tagDOF(i),'r') %Revolute joint
%         %Linear
%         if i == 1
%             Jv(:,i) = cross(ziminus1_w(:,i),ti_w(:,numDOF));
%         else
%             Jv(:,i) = cross(ziminus1_w(:,i),ti_w(:,numDOF)-ti_w(:,i-1));
%         end
%         %Angular
%         Jw(:,i) = ziminus1_w(:,i);
%     elseif strcmp(tagDOF(i),'p') %Prismatic joint
%         Jv(:,i) = ziminus1_w(:,i);
%     else
%         disp("Joint Tags not Proper");
%     end
% end
% 
% 
% Jn_w = [Jv;Jw];
% Jn_w = simplify(Jn_w);
% 
% %% Calculate Jacobians for Center of Masses
% 
% Jcmi_w = sym(zeros(6,numDOF,numDOF));
% tcmi_w = sym(zeros(3,numDOF));
% for i = 1:numDOF
%     tcmi_w(:,i) =   Hcmi_w(1:3,4,i);
% end
% 
% for i = 1:numDOF
% 
% Jv = sym(zeros(3,numDOF));
% Jw = sym(zeros(3,numDOF));
%     for j = 1:i
%         if strcmp(tagDOF(j),'r') %Revolute Joint
%             %Linear
%             if j == 1
%             Jv(:,j) = cross(ziminus1_w(:,j),tcmi_w(:,i));
%             else
%             Jv(:,j) = cross(ziminus1_w(:,j),tcmi_w(:,i)-ti_w(:,j-1));    
%             end
%             %Angular
%             Jw(:,j) = ziminus1_w(:,j);
%             
%         elseif strcmp(tagDOF(j),'p') %Prismatic Joint
%             Jv(:,j) = ziminus1_w(:,j);
%    
%         else
%             disp("Joint Tags not Proper");
%         end
%     end     
%    
% Jcmi_w(:,:,i) = simplify([Jv;Jw]);
% end




%%  Define axis and radius wrt. base
% kz = [0;0;1];
% ziminus1_0 = sym(zeros(3,numDOF));
% ti_0 = sym(zeros(3,numDOF));
% 
% ziminus1_0(:,1) = kz;
% for i = 2:numDOF
%    ziminus1_0(:,i) = Hi_0(1:3,1:3,i-1)*kz;
% end
% 
% for i = 1:numDOF
%    ti_0(:,i) = Hi_0(1:3,4,i);
% end


%% Calculate End Effector Jacobian

% Jv = sym(zeros(3,numDOF));
% Jw = sym(zeros(3,numDOF));
% 
% 
% for i = 1:numDOF
%     if strcmp(tagDOF(i),'r') %Revolute joint
%         %Linear
%         if i == 1
%             Jv(:,i) = cross(ziminus1_0(:,i),ti_0(:,numDOF));
%         else
%             Jv(:,i) = cross(ziminus1_0(:,i),ti_0(:,numDOF)-ti_0(:,i-1));
%         end
%         %Angular
%         Jw(:,i) = ziminus1_0(:,i);
%     elseif strcmp(tagDOF(i),'p') %Prismatic joint
%         Jv(:,i) = ziminus1_0(:,i);
%     else
%         disp("Joint Tags not Proper");
%     end
% end
% 
% Jn_0 = [Jv;Jw];
% Jn_0 = simplify(Jn_0);
% 
% 
% dJn_0 = sym(zeros(6,numDOF,numDOF));
% for i = 1:numDOF
%     dJn_0(:,:,i)=  simplify(diff(Jn_0(:,:),q(i)));
%     
% end
% 
% pJn_0 = sym(zeros(6,numDOF));
% 
% for i = 1:numDOF
%     pJn_0=  simplify(pJn_0 +dJn_0(:,:,i)*qip(i));
%     
% end


%% Calculate Jacobians for Center of Masses

% Jcmi_0 = sym(zeros(6,numDOF,numDOF));
% tcmi_0 = sym(zeros(3,numDOF));
% for i = 1:numDOF
%     tcmi_0(:,i) =   Hcmi_0(1:3,4,i);
% end
% 
% for i = 1:numDOF
% 
% Jv = sym(zeros(3,numDOF));
% Jw = sym(zeros(3,numDOF));
%     for j = 1:i
%         if strcmp(tagDOF(j),'r') %Revolute Joint
%             %Linear
%             if j == 1
%             Jv(:,j) = cross(ziminus1_0(:,j),tcmi_0(:,i));
%             else
%             Jv(:,j) = cross(ziminus1_0(:,j),tcmi_0(:,i)-ti_0(:,j-1));    
%             end
%             %Angular
%             Jw(:,j) = ziminus1_0(:,j);
%             
%         elseif strcmp(tagDOF(j),'p') %Prismatic Joint
%             Jv(:,j) = ziminus1_0(:,j);
%    
%         else
%             disp("Joint Tags not Proper");
%         end
%     end     
%    
% Jcmi_0(:,:,i) = simplify([Jv;Jw]);
% end
% 
% %% Write Jacobians to file
% 
% 
% fid = fopen('MyJacobmatrixUR10ProjectPose.txt','wt');
% WriteToFile(fid,Jn_w,'Jn_w');
% for i = 1:numDOF
%    WriteToFile(fid,Jcmi_w(:,:,i),strcat('Jcm',num2str(i),'_w')) 
% end
% 
% WriteToFile(fid,Jn_0,'Jn_0');
% for i = 1:numDOF
%    WriteToFile(fid,Jcmi_0(:,:,i),strcat('Jcm',num2str(i),'_0')) 
% end
% 
% WriteToFile(fid,dJn_0,'dJn_0');
% 
% WriteToFile(fid,pJn_0,'pJn_0');
% 
% fclose(fid);


%% Dynamics

syms I111 I112 I113 I122 I123 I133 I211 I212 I213 I222 I223 I233 I311 I312 I313 I322 I323 I333 real
syms I411 I412 I413 I422 I423 I433 I511 I512 I513 I522 I523 I533 I611 I612 I613 I622 I623 I633 real
syms gx gy gz g real
syms m1 m2 m3 m4 m5 m6 real

%%

Inertias = [I111 I112 I113 I122 I123 I133 I211 I212 I213 I222 I223 I233 I311 I312 I313 I322 I323 I333];
Inertias = [Inertias,[I411 I412 I413 I422 I423 I433 I511 I512 I513 I522 I523 I533 I611 I612 I613 I622 I623 I633]];
I = sym(zeros(3,3,numDOF));

I(:,:,1) = [I111,I112,I113;
    I112,I122,I123;
    I113,I123,I133];
I(:,:,2) = [I211,I212,I213;
    I212,I222,I223;
    I213,I223,I233];
I(:,:,3) = [I311,I312,I313;
    I312,I322,I323;
    I313,I323,I333];
I(:,:,4) = [I411,I412,I413;
    I412,I422,I423;
    I413,I423,I433];
I(:,:,5) = [I511,I512,I513;
    I512,I522,I523;
    I513,I523,I533];
I(:,:,6) = [I611,I612,I613;
    I612,I622,I623;
    I613,I623,I633];
m = [m1,m2,m3,m4,m5,m6];
gs = [gx,gy,gz]*g;



M = sym(zeros(numDOF,numDOF));
for i = 1:numDOF
    M = M + m(i)*Jcmi_0(1:3,:,i)'*Jcmi_0(1:3,:,i) +...
        Jcmi_0(4:6,:,i)'*Hcmi_0(1:3,1:3,i)*I(:,:,i)*Hcmi_0(1:3,1:3,i)'*Jcmi_0(4:6,:,i);
end
M = simplify(M);

% M = sym(zeros(3,3));
% for i = 1:3
%     M = M + m(i)*Jcmi_0(1:3,:,i)'*Jcmi_0(1:3,:,i) +...
%         Jcmi_0(4:6,:,i)'*Hcmi_0(1:3,1:3,i)*I(:,:,i)*Hcmi_0(1:3,1:3,i)'*Jcmi_0(4:6,:,i);
% end
% M = simplify(M);

dMqi = sym(zeros(numDOF,numDOF,numDOF));
for i = 1:numDOF
    dMqi(:,:,i)= simplify(diff(M,q(i)));
    
end


C = sym(zeros(numDOF,numDOF));
for k = 1:numDOF
    for j = 1:numDOF
        indexC = sym(zeros(1));
        for i = 1:numDOF
            indexC = indexC + (1/2)*(dMqi(k,j,i) +dMqi(k,i,j) - dMqi(i,j,k))*qip(i);
        end
        indexC = simplify(indexC);
        C(k,j) = indexC;
    end
end
C = simplify(C);

P = sym(zeros(1));
for i = 1:numDOF
   P = P + m(i)*(gs*tcmi_0(:,i));
end
P = simplify(P);

G = sym(zeros(numDOF,1));
for i = 1:numDOF
    G(i,1) = simplify(diff(P,q(i)));
end

Mp = sym(zeros(numDOF,numDOF));

for i = 1:numDOF
   Mp = Mp + dMqi(:,:,i)*qip(i); 
end
Mp = simplify(Mp);
% 
% ValidateMatProp(M,Mp,C,[L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,m1,m2,m3,...
%     I111 I112 I113 I122 I123 I133 I211 I212 I213 I222 I223 I233 I311 I312 I313 I322 I323 I333],...
%     [0.128,0.1639,0.6127,0.0922,0.1,0.15,0.30635,0.1,0.34365,[5.0;15.0;15.0]',...
%     [0.04;0.02;0.02;0.04;0.02;0.04]',[0.04;0.02;0.02;0.04;0.02;0.04]',...
%     [0.04;0.02;0.02;0.04;0.02;0.04]'],[q1,q2,q3]);

ValidateMatProp(M,Mp,C,[L, m,Inertias],...
        [0.128,0.1639,0.6127,0.0922,0.1,0.15,0.30635,0.1,0.34365,[5.0;15.0;15.0]',...
    [0.04;0.02;0.02;0.04;0.02;0.04]',[0.04;0.02;0.02;0.04;0.02;0.04]',...
    [0.04;0.02;0.02;0.04;0.02;0.04]'],q,numDOF);

% MyM = M;
% MyC = C;
% MyG = G;
% MyMp = Mp;
%%
syms q1pp q2pp q3pp q4pp q5pp q6pp real

%%
qipp = [q1pp,q2pp,q3pp,q4pp,q5pp,q6pp];
Tau = expand(M*qipp' + C*qip' + G);

fid = fopen('MyDynEqnUR10Project.txt','wt');
for i = 1:numDOF
    for j = 1:numDOF
        WriteToFile(fid,M(i,j),strcat('M(',num2str(i),',',num2str(j),')'),'noSpace');        
    end
end
fprintf(fid,'\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid,'\n');
fprintf(fid,'\n');
for i = 1:numDOF
    for j = 1:numDOF
        WriteToFile(fid,C(i,j),strcat('C(',num2str(i),',',num2str(j),')'),'noSpace');
    end
end
fprintf(fid,'\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid,'\n');
fprintf(fid,'\n');
for i = 1:numDOF
    WriteToFile(fid,G(i),strcat('G(',num2str(i),',1)'),'noSpace');
end
fprintf(fid,'\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid,'\n');
fprintf(fid,'\n');
for i = 1:numDOF
    for j = 1:numDOF
        WriteToFile(fid,Mp(i,j),strcat('Mp(',num2str(i),',',num2str(j),')'),'noSpace');
    end
end
fprintf(fid,'\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid,'\n');
fprintf(fid,'\n');
for i = 1:numDOF
    WriteToFile(fid,Tau(i),strcat('Tau(',num2str(i),',1)'),'noSpace');
end
% WriteToFile(fid,M,'M');
% WriteToFile(fid,C,'C');
% WriteToFile(fid,G,'G');
% WriteToFile(fid,Mp,'Mp');
% WriteToFile(fid,Tau,'Tau');
fclose(fid);


%% Calculate the regressor model

% Regressor with all parameters split
vars_all = cell(numDOF,1);
coef_all = cell(numDOF,1);
for i = 1:numDOF 
   [coef_all{i}, vars_all{i}] = coeffs(Tau(i),[m L Inertias]); 
end


unique_vars_all = unique([vars_all{1},vars_all{2},vars_all{3},vars_all{4},vars_all{5},vars_all{6}]);

Theta = unique_vars_all';
var_num_all = size(Theta,1);
Y = sym(zeros(numDOF,var_num_all));

for i = 1:var_num_all
   for j = 1:numDOF
       sz = size(vars_all{j},2);
        for k = 1:sz
            if isequaln(vars_all{j}(k),unique_vars_all(i))
                Y(j,i) = coef_all{j}(k);
            end
        end
   end
end

assert(isequaln(expand(Y*Theta),Tau),'Regressor not correct');
% TODO
% Regressor with compact parameters
% vars_comp = cell(3,1);
% coef_comp = cell(3,1);
% [coef_comp{1}, vars_comp{1}] = coeffs(Tau(1),[q1 q2 q3 q1p q2p q3p q1pp q2pp q3pp cos(q1) sin(q1) cos(q2) sin(q2) cos(q3) sin(q3)]);
% [coef_comp{2}, vars_comp{2}] = coeffs(Tau(2),[q1 q2 q3 q1p q2p q3p q1pp q2pp q3pp cos(q1) sin(q1) cos(q2) sin(q2) cos(q3) sin(q3)]);
% [coef_comp{3}, vars_comp{3}] = coeffs(Tau(3),[q1 q2 q3 q1p q2p q3p q1pp q2pp q3pp cos(q1) sin(q1) cos(q2) sin(q2) cos(q3) sin(q3)]);


%% Calculate the reference regressor model

syms q1pr q2pr q3pr q4pr q5pr q6pr real
syms q1ppr q2ppr q3ppr q4ppr q5ppr q6ppr real
%%
qipr = [q1pr q2pr q3pr q4pr q5pr q6pr];
qippr = [q1ppr q2ppr q3ppr q4ppr q5ppr q6ppr];

Taur = expand(M*qippr' + C*qipr' + G);



% Regressor with all parameters split
vars_all = cell(numDOF,1);
coef_all = cell(numDOF,1);
% [coef_all{1}, vars_all{1}] = coeffs(Taur(1),[ m1 m2 m3 L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 I111 I112 I113 I122 I123 I133 I211 I212 I213 I222 I223 I233 I311 I312 I313 I322 I323 I333]);
% [coef_all{2}, vars_all{2}] = coeffs(Taur(2),[ m1 m2 m3 L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 I111 I112 I113 I122 I123 I133 I211 I212 I213 I222 I223 I233 I311 I312 I313 I322 I323 I333]);
% [coef_all{3}, vars_all{3}] = coeffs(Taur(3),[ m1 m2 m3 L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 I111 I112 I113 I122 I123 I133 I211 I212 I213 I222 I223 I233 I311 I312 I313 I322 I323 I333]);
for i = 1:numDOF
    [coef_all{i}, vars_all{i}] = coeffs(Taur(i),[m L Inertias]);
end

unique_vars_all = unique([vars_all{1},vars_all{2},vars_all{3},vars_all{4},vars_all{5},vars_all{6}]);

Thetar = unique_vars_all';
var_num_all = size(Thetar,1);
Yr = sym(zeros(numDOF,var_num_all));

for i = 1:var_num_all
   for j = 1:numDOF
       sz = size(vars_all{j},2);
        for k = 1:sz
            if isequaln(vars_all{j}(k),unique_vars_all(i))
                Yr(j,i) = coef_all{j}(k);
            end
        end
   end
end

assert(isequaln(expand(Yr*Thetar),Taur),'Regressor not correct');


fid = fopen('MyRegressorUR10Project.txt','wt');
for i = 1:numDOF
    WriteToFile(fid,Tau(i,1),strcat('Tau(',num2str(i),',1)'),'noSpace');
end
fprintf(fid,'\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid,'\n');
fprintf(fid,'\n');
for i = 1:numDOF
    for j = 1:var_num_all
        WriteToFile(fid,Y(i,j),strcat('Y(',num2str(i),',',num2str(j),')'),'noSpace');
    end
end
fprintf(fid,'\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid,'\n');
fprintf(fid,'\n');
% WriteToFile(fid,Theta,'Theta');

% WriteToFile(fid,Y,'Y','verbose');
WriteToFile(fid,Theta,'Theta','verbose');

fclose(fid);

fid = fopen('MyRegressorUR10refProject.txt','wt');
for i = 1:numDOF
    WriteToFile(fid,Taur(i,1),strcat('Taur(',num2str(i),',1)'),'noSpace');
end
fprintf(fid,'\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid,'\n');
fprintf(fid,'\n');
for i = 1:numDOF
    for j = 1:var_num_all
        WriteToFile(fid,Yr(i,j),strcat('Yr(',num2str(i),',',num2str(j),')'),'noSpace');
    end
end
fprintf(fid,'\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
fprintf(fid,'\n');
fprintf(fid,'\n');
% WriteToFile(fid,Thetar,'Theta');

% WriteToFile(fid,Yr,'Yr','verbose');
WriteToFile(fid,Thetar,'Theta','verbose');

fclose(fid);
