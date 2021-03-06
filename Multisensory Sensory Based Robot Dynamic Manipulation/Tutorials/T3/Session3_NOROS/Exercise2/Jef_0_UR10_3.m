function Jef = Jef_0_UR10_3(in1)
%JEF_0_UR10_3
%    JEF = JEF_0_UR10_3(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    21-Feb-2018 21:13:19

% J(q)
q1 = in1(1,:);
q2 = in1(2,:);
q3 = in1(3,:);
Jef = reshape([cos(q1).*2.561e-1+sin(q1).*sin(q2).*6.127e-1-cos(q2).*sin(q1).*sin(q3).*6.873e-1+cos(q3).*sin(q1).*sin(q2).*6.873e-1,sin(q1).*2.561e-1-cos(q1).*sin(q2).*6.127e-1+cos(q1).*cos(q2).*sin(q3).*6.873e-1-cos(q1).*cos(q3).*sin(q2).*6.873e-1,0.0,0.0,0.0,1.0,-cos(q1).*(cos(q2-q3).*6.873e-1+cos(q2).*6.127e-1),-sin(q1).*(cos(q2-q3).*6.873e-1+cos(q2).*6.127e-1),sin(q2-q3).*(-6.873e-1)-sin(q2).*6.127e-1,sin(q1),-cos(q1),0.0,cos(q2-q3).*cos(q1).*6.873e-1,cos(q2-q3).*sin(q1).*6.873e-1,sin(q2-q3).*6.873e-1,-sin(q1),cos(q1),0.0],[6,3]);
