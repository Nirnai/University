function [ tau ] = Tau( u )
%TAU Summary of this function goes here
%   Detailed explanation goes here
% Times are defined in TrajGen.m file
global t0 t1 t2 t3 t4 t5 t6
persistent Q_error_total
t=u(1);

q1d=deg2rad(u(2));
q2d=deg2rad(u(3));
q3d=deg2rad(u(4));

q1dp=deg2rad(u(5));
q2dp=deg2rad(u(6));
q3dp=deg2rad(u(7));

q1=u(14);
q2=u(15);
q3=u(16);

q1p=u(17);
q2p=u(18);
q3p=u(19);

Q=[q1;q2;q3];
Qp=[q1p;q2p;q3p];

Qd=[q1d;q2d;q3d];
Qdp=[q1dp;q2dp;q3dp];

%Joint Errors
Q_error = Q -Qd;
Qp_error = Qp-Qdp;
%Robot Parameters

m1=u(20);
m2=u(21);
m3=u(22);
g=u(23);

L1=u(24);
L2=u(25);
L4=u(26);
L6=u(27);
L7=u(28);
L9=u(29);
L3=u(30);
L5=u(31);
L8=u(32);
L10=u(33);






Kd=diag([u(8);u(9);u(10)]);
Kp=diag([u(11);u(12);u(13)]);
Ki=diag([u(37);u(38);u(39)]);




gx=u(34);
gy=u(35);
gz=u(36);


G = [...
- g*gy*m3*(sin(q1)*(L9 + 3/200) - L3*cos(q1)*sin(q2) + L10*cos(q1)*cos(q2)*sin(q3) - L10*cos(q1)*cos(q3)*sin(q2)) - g*gx*m3*(cos(q1)*(L9 + 3/200) + L3*sin(q1)*sin(q2) - L10*cos(q2)*sin(q1)*sin(q3) + L10*cos(q3)*sin(q1)*sin(q2)) - g*gx*m2*(L7*cos(q1) + L8*sin(q1)*sin(q2)) - g*gy*m2*(L7*sin(q1) - L8*cos(q1)*sin(q2));...
            g*gx*m3*(L3*cos(q1)*cos(q2) + L10*cos(q1)*cos(q2)*cos(q3) + L10*cos(q1)*sin(q2)*sin(q3)) - g*gz*m3*(L3*sin(q2) + L10*sin(q2 - q3)) + g*gy*m3*(L3*cos(q2)*sin(q1) + L10*cos(q2)*cos(q3)*sin(q1) + L10*sin(q1)*sin(q2)*sin(q3)) - L8*g*gz*m2*sin(q2) + L8*g*gx*m2*cos(q1)*cos(q2) + L8*g*gy*m2*cos(q2)*sin(q1);...
                                                                                                                                                    L10*g*gz*m3*sin(q2 - q3) - g*gy*m3*(L10*cos(q2)*cos(q3)*sin(q1) + L10*sin(q1)*sin(q2)*sin(q3)) - g*gx*m3*(L10*cos(q1)*cos(q2)*cos(q3) + L10*cos(q1)*sin(q2)*sin(q3))];
 

%Define controllers PD, PD+G, PID
if(t <= t0)
    % PD
    tauc = -Kp*Q_error - Kd*Qp;
elseif(t>=t0 && t<t1)
    % PD+G
    tauc = -Kp*Q_error - Kd*Qp + G;
elseif(t>=t1 && t<t2)
    % PID
    if(t == t1)
        Q_error_total = [0;0;0];
    end    

    Q_error_total = Q_error_total + Q_error;
    tauc = -Kp*Q_error - Kd*Qp_error - Ki * Q_error_total;
elseif(t>=t2 && t<t3)
    % PD
    tauc = -Kp*Q_error - Kd*Qp;
elseif(t>=t3 && t<t4)
    % PD+G
    tauc = -Kp*Q_error - Kd*Qp + G;
else
    % PID
    if(t == t4)
        Q_error_total = [0;0;0];
    end    

    Q_error_total = Q_error_total + Q_error;
    tauc = -Kp*Q_error - Kd*Qp_error - Ki * Q_error_total;
end
DeltaQ = Q_error;

DeltaQp = Qp_error;

tau=[tauc;DeltaQ;DeltaQp];

end

