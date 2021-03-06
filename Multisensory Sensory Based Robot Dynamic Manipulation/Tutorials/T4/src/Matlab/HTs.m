function [H1_W,H2_W,H3_W] = HTs(in1)
%HTS
%    [H1_W,H2_W,H3_W] = HTS(IN1)

%    This function was generated by the Symbolic Math Toolbox version 8.0.
%    30-Dec-2017 12:53:03

l1 = in1(:,7);
l2 = in1(:,8);
l3 = in1(:,9);
q1 = in1(:,1);
q2 = in1(:,2);
q3 = in1(:,3);
H1_W = reshape([1.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,-1.0,0.0,0.0,0.0,0.0,-l1-q1,1.0],[4,4]);
if nargout > 1
    t2 = sin(q2);
    t3 = cos(q2);
    t4 = q2+q3;
    t5 = sin(t4);
    t6 = l2.*t2;
    H2_W = reshape([-t2,0.0,t3,0.0,-t3,0.0,-t2,0.0,0.0,-1.0,0.0,0.0,t6,0.0,-l1-q1-l2.*t3,1.0],[4,4]);
end
if nargout > 2
    t7 = cos(t4);
    H3_W = reshape([-t5,0.0,t7,0.0,-t7,0.0,-t5,0.0,0.0,-1.0,0.0,0.0,t6+l3.*t5,0.0,-l1-q1-l2.*t3-l3.*t7,1.0],[4,4]);
end
