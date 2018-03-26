function  [HT_Joints, HT_CM] = DH2H(DH_Joints, DH_CenterOfMasses, u)

% DH2H  Converts DH-Table to Homogenous Tansformations.
%       The function generates functions which can be called with
%       numeric inputs u.
%
%   H = DH2H(DH,u)          generates relative Homogenous transformations 
%                           between joints.
%   H = DH2H(DH,u,true)     generates relative Homogenous transformations 
%                           between the previous joints and center of masses.

%%


[n,~] = size(DH_Joints);
[m,~] = size(DH_CenterOfMasses);
HT_Joints = cell(2,n);
HT_CM = cell(2,m);


% Relative Transformations Joints
for i = 1:n
    D = DH_Joints(i,:);
    H = [...
     cos(D(1)), -sin(D(1)) * cos(D(4)),  sin(D(1)) * sin(D(4)), D(3) * cos(D(1)); ...
     sin(D(1)),  cos(D(1)) * cos(D(4)), -cos(D(1)) * sin(D(4)), D(3) * sin(D(1)); ...
             0,              sin(D(4)),              cos(D(4)),             D(2); ...
             0,                      0,                      0,               1];
   HT_Joints{1,i} = simplify(H);
end

% Relative Transformations Center of Masses

for i = 1:m
    D = DH_CenterOfMasses(i,:);
    H = [...
     cos(D(1)), -sin(D(1)) * cos(D(4)),  sin(D(1)) * sin(D(4)), D(3) * cos(D(1)); ...
     sin(D(1)),  cos(D(1)) * cos(D(4)), -cos(D(1)) * sin(D(4)), D(3) * sin(D(1)); ...
             0,              sin(D(4)),              cos(D(4)),             D(2); ...
             0,                      0,                      0,               1];
   HT_CM{1,i} = simplify(H);
end


% Absolute Transformations to Base
HT_Joints{2,1} = HT_Joints{1,1};
HT_CM{2,1} = HT_CM{1,1};
for i = 2:n
    H = HT_Joints{2,1};
    for j = 2:i
       H = H * HT_Joints{1,j};
    end
    HT_Joints{2,i} = simplify(H);
    if(i<=m)
        Hcm = HT_Joints{2,i-1} * HT_CM{1,i};
        HT_CM{2,i} = simplify(Hcm);
    end
end





