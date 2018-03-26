function  [HT] = DH2H(DH)

% DH2H  Converts DH-Table to Homogenous Tansformations.
%       The function generates functions which can be called with
%       numeric inputs u.
%
%   H = DH2H(DH,u)          generates relative Homogenous transformations 
%                           between joints.
%   H = DH2H(DH,u,true)     generates relative Homogenous transformations 
%                           between the previous joints and center of masses.

%%


[n,~] = size(DH);
HT = cell(2,n);

% Relative Transformations Joints
for i = 1:n
    D = DH(i,:);
    H = [...
     cos(D(1)), -sin(D(1)) * cos(D(4)),  sin(D(1)) * sin(D(4)), D(3) * cos(D(1)); ...
     sin(D(1)),  cos(D(1)) * cos(D(4)), -cos(D(1)) * sin(D(4)), D(3) * sin(D(1)); ...
             0,              sin(D(4)),              cos(D(4)),             D(2); ...
             0,                      0,                      0,               1];
         
   HT{1,i} = simplify(H);
end


%% Absolute Transformations to Base
HT{2,1} = HT{1,1};
HT{2,2} = HT{2,1} * HT{1,2};
HT{2,3} = HT{2,2} * HT{1,3};





