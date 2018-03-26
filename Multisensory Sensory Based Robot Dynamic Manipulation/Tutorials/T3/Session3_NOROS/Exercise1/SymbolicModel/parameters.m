function t = parameters(t)

global L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11 L12 L13 L14 L15 L16 L17 L18 c T1 T3 

L1  = 0.300; 
L2  = 0.340; 
L3  = 0.583; 
L4  = 0.080; 
L5  = 0.080 + 0.15;  % Including the length of Joint
L6  = 0.115;
L7  = 0.082 + 0.15;  % Including the length of Joint 
L8  = 0.108; 
L9  = 0.123; 
L10 = 0.167; 
L11 = 0.028;
L12 = 0.050; 

c   = 0.320;
T1  =  pi/4; 

L13 = L9 * sin(T1);
L14 = L9 * cos(T1);
L15 = (c/2 + L12) * cos(T1);
L16 = (c/2 + L12) * sin(T1);
L17 = L10 * 1/cos(T1);

T3 = atan(L14/L17);

L18 = L17 * 1/cos(T3);


end

