function V = Lyapunov(u)
persistent l
t=u(39);
%Joint Position
q1=u(1);
q2=u(2);
q3=u(3);

Q = [q1;q2;q3];
%Joint Velocity
qp1=u(4);
qp2=u(5);
qp3=u(6);

Qp = [qp1;qp2;qp3];
%Kinematic Parameters

L7=u(11);
L9=u(12);
L3=u(13);
L8=u(15);
L10=u(16);

%Dynamic Parameters
m2=u(18);
m3=u(19);


I133=u(25);

I211=u(26);
I212=u(27);
I213=u(28);
I222=u(29);
I223=u(30);
I233=u(31);

I311=u(32);
I312=u(33);
I313=u(34);
I322=u(35);
I323=u(36);
I333=u(37);



q1d=deg2rad(u(56));
q2d=deg2rad(u(57));
q3d=deg2rad(u(58));

Qd = [q1d;q2d;q3d];

Q_error = Q - Qd;
%Inertia Matrix
M = [... 
  I133 + I211/2 + I222/2 + I311/2 + I322/2 + (9*m3)/40000 + (I311*cos(2*q2 - 2*q3))/2 - (I322*cos(2*q2 - 2*q3))/2 + (3*L9*m3)/100 + I312*sin(2*q2 - 2*q3) + (L3^2*m3)/2 + L7^2*m2 + (L8^2*m2)/2 + L9^2*m3 + (L10^2*m3)/2 + (I211*cos(2*q2))/2 - (I222*cos(2*q2))/2 - I212*sin(2*q2) - (L3^2*m3*cos(2*q2))/2 - (L8^2*m2*cos(2*q2))/2 - (L10^2*m3*cos(2*q2 - 2*q3))/2 + L3*L10*m3*cos(q3) - L3*L10*m3*cos(2*q2 - q3), I213*cos(q2) - I223*sin(q2) - I313*cos(q2 - q3) - I323*sin(q2 - q3) - (3*L10*m3*cos(q2 - q3))/200 - (3*L3*m3*cos(q2))/200 - L3*L9*m3*cos(q2) - L7*L8*m2*cos(q2) - L9*L10*m3*cos(q2 - q3), I313*cos(q2 - q3) + I323*sin(q2 - q3) + (3*L10*m3*cos(q2 - q3))/200 + L9*L10*m3*cos(q2 - q3);...
                                                                                                                                                                                                                          I213*cos(q2) - I223*sin(q2) - I313*cos(q2 - q3) - I323*sin(q2 - q3) - (3*L10*m3*cos(q2 - q3))/200 - (3*L3*m3*cos(q2))/200 - L3*L9*m3*cos(q2) - L7*L8*m2*cos(q2) - L9*L10*m3*cos(q2 - q3),                                                                                                                         m3*L3^2 + 2*m3*cos(q3)*L3*L10 + m2*L8^2 + m3*L10^2 + I233 + I333,                                                        - m3*L10^2 - L3*m3*cos(q3)*L10 - I333;...
                                                                                                                                                                                                                                                                                                                      I313*cos(q2 - q3) + I323*sin(q2 - q3) + (3*L10*m3*cos(q2 - q3))/200 + L9*L10*m3*cos(q2 - q3),                                                                                                                                                    - m3*L10^2 - L3*m3*cos(q3)*L10 - I333,                                                                              m3*L10^2 + I333];
Kp=diag([u(65);u(66);u(67)]);



%% Lyapunov Function
V = 1/2 * Qp' * M * Qp + 1/2 * Q_error' * Kp * Q_error;


%% Phasendiagramm
% if(t == 0)
%     figure(1);
%     axis([-1 1 -1 1 -1 1]) 
%     l = animatedline;
% end
% grid on;
% addpoints(l,q1-q1d,q2-q2d,q3-q3d);
% drawnow limitrate
end

