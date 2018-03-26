function [qd] = TrajGen_JointSpace( u )
% %TRAJGEN Summary of this function goes here
% %   Detailed explanation goes here
% 
global t0 t1 t2 t3


t0 = 10;
t1 = 20;
t2 = 40;
t3 = 60;

PID = false;

Ki = diag([u(7);u(8);u(9)]);

if(~isempty(nonzeros(Ki)))
    PID = true;
end


q1 = u(4);
q2 = u(5);
q3 = u(6);
Q = [q1;q2;q3];

t = u(10);


tref = 1; 

A = 0.5;
w1 = 1;
w2 = 2;

persistent Q_init t_old tr

dt = 0;

if (t == 0)
    t_old = t;
else
    dt = t - t_old;
    t_old = t;
end
    
if(t<t0)   
    if(t == 0)
        Q_init = Q;
        tr = 0;
    end
    tr = tr + dt;
    
    qd1 = deg2rad(u(1));
    qd2 = deg2rad(u(2));
    qd3 = deg2rad(u(3));
    qdp1 = 0;
    qdp2 = 0;
    qdp3 = 0;
    qdpp1 = 0;
    qdpp2 = 0;
    qdpp3 = 0;
elseif(t>=t0 && t<t1)
    if(t == t0)
        Q_init = Q;
        tr = 0;
    end
    tr = tr + dt;
    
    qd1 = 0;
    qd2 = 0;
    qd3 = 0;
    qdp1 = 0;
    qdp2 = 0;
    qdp3 = 0;
    qdpp1 = 0;
    qdpp2 = 0;
    qdpp3 = 0;
elseif(t>=t1 && t<t2)
    if(t == t1)
        Q_init = Q;
        tr = 0;
    end
    
    tr = tr + dt;
    
    qd1 = A*sin(w1*t);
    qd2 = A*sin(w1*t);
    qd3 = A*sin(w1*t);
    qdp1 = A*w1*cos(w1*t);
    qdp2 = A*w1*cos(w1*t);
    qdp3 = A*w1*cos(w1*t);
    qdpp1 = -A*w1^2*sin(w1*t);
    qdpp2 = -A*w1^2*sin(w1*t);
    qdpp3 = -A*w1^2*sin(w1*t);
else
    
    if(t == t2)
        Q_init = Q;
        tr = 0;
    end
    
    tr = tr + dt;
    
    qd1 = A*sin(w2*t);
    qd2 = A*sin(w2*t);
    qd3 = A*sin(w2*t);
    qdp1 = A*w2*cos(w2*t);
    qdp2 = A*w2*cos(w2*t);
    qdp3 = A*w2*cos(w2*t);
    qdpp1 = -A*w2^2*sin(w2*t);
    qdpp2 = -A*w2^2*sin(w2*t);
    qdpp3 = -A*w2^2*sin(w2*t);
end

%% Reduce Windup in PID Controller


if(PID)
    if(abs(qd1 - Q_init(1)) > 0.01)
        qd1 = qd1 - exp(-tr/tref)*(qd1-Q_init(1));
        qdp1 = -1/tref * exp(-tr/tref)*(qd1-Q_init(1));
        qdpp1 = 1/tref^2 * exp(-tr/tref)*(qd1-Q_init(1));
    end

    if(abs(qd2 - Q_init(2)) > 0.01)
        qd2 = qd2 - exp(-tr/tref)*(qd2-Q_init(2));
        qdp2 = -1/tref * exp(-tr/tref)*(qd2-Q_init(2));
        qdpp2 = -1/tref^2 * exp(-tr/tref)*(qd2-Q_init(2));
    end

    if(abs(qd3 - Q_init(3)) > 0.01)
        qd3 = qd3 - exp(-tr/tref)*(qd3-Q_init(3));
        qdp3 = -1/tref * exp(-tr/tref)*(qd3-Q_init(3));
        qdpp3 = -1/tref^2 * exp(-tr/tref)*(qd3-Q_init(3));
    end
end
% 
% 
qd=[qd1;qd2;qd3;qdp1;qdp2;qdp3;qdpp1;qdpp2;qdpp3];
% qd=zeros(9,1);

end

