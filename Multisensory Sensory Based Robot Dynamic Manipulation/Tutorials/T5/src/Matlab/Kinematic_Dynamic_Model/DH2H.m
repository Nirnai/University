function  [HT] = DH2H(DH)

% DH2H  Converts DH-Table to Homogenous Tansformations.
%       The function generates functions which can be called with
%       numeric inputs u.
%
%   H = DH2H(DH,u)          generates relative Homogenous transformations 
%                           between joints.
%   H = DH2H(DH,u,true)     generates relative Homogenous transformations 
%                           between the previous joints and center of masses.


[n,~] = size(DH);

% Initialize HT
HT = cell(4,n/2);

% HT holds the transformations in following order:
% HT = { H1_0     H2_1   H3_2 ... HN_N-1 
%        H1_0     H2_0   H3_0 ... HN_0
%        Hcm1_0 Hcm2_1 Hcm3_2 ... HcmN_N-1
%        Hcm1_0 Hcm2_0 Hcm3_0 ... HcmN_0 }



% Relative Transformations Joints
for i = 1:n/2
    % Compute Transformation to Joint from DH-Row
    JT = DH(i,:);
    CM = DH(i+n/2,:);
    D = [JT;CM];
    for j = 1:size(D,1)
        H = [...
         cos(D(j,1)), -sin(D(j,1)) * cos(D(j,4)),  sin(D(j,1)) * sin(D(j,4)), D(j,3) * cos(D(j,1)); ...
         sin(D(j,1)),  cos(D(j,1)) * cos(D(j,4)), -cos(D(j,1)) * sin(D(j,4)), D(j,3) * sin(D(j,1)); ...
                   0,                sin(D(j,4)),                cos(D(j,4)),               D(j,2); ...
                   0,                          0,                          0,                   1];
               
       if(j == 1)        
        HT{1,i} = simplify(H);
       else
        HT{3,i} = simplify(H);
       end
   end
end


% Absolute Transformations to Base
HT{2,1} = HT{1,1};
HT{2,2} = simplify(HT{2,1} * HT{1,2});
HT{2,3} = simplify(HT{2,2} * HT{1,3});

HT{4,1} = HT{3,1};
HT{4,2} = simplify(HT{2,1} * HT{3,2});
HT{4,3} = simplify(HT{2,2} * HT{3,3});




