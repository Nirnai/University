function [Out] = TrajGen_OperationalSpace( u )
t = u(10);

q1 = u(4);
q2 = u(5);
q3 = u(6);

L1 = u(11);
L2 = u(12);
L3 = u(13);
L4 = u(14);
L5 = u(15);

Q_Current = zeros(9,1);
Q_Current(1:3) = [q1;q2;q3];

X_Current = zeros(9,1);
X_Current(1:3) = [sin(q1)*(L2 + L4) - L3*cos(q1)*cos(q2) - L5*cos(q1)*cos(q2)*cos(q3) + L5*cos(q1)*sin(q2)*sin(q3); 
     L5*sin(q1)*sin(q2)*sin(q3) - L3*cos(q2)*sin(q1) - L5*cos(q2)*cos(q3)*sin(q1) - cos(q1)*(L2 + L4); 
     L1 - L5*sin(q2 + q3) - L3*sin(q2)];

persistent a b c

% States
if(t == 0)
    a = true; 
    b = false; 
    c = false;
end

Q_Home = zeros(9,1);
Q_Home(2) = -pi/2;
Q_Home(3) =  pi/2;

if(a)
   Q_Desired = Q_Home;
   X_Desired = zeros(9,1);
   if( norm(Q_Current - Q_Desired)< 1e-2 )
       a = false;
       b = true;
   end
end

if(b)
    Q_Desired = Q_Home;
    X_Desired = zeros(9,1);
    X_Desired(1) = -0.75;
    X_Desired(2) = 0.5;
    X_Desired(3) = 0.1;
    if( norm(X_Current - X_Desired)< 0.4 )
        a = false;
        b = false;
        c = true;
    end
end

if(c)
     Q_Desired = Q_Home;
     X_Desired =  zeros(9,1);
% 
%     w1 = 1;
%     w2 = 5;
%     
%     X_Desired(1) = 0.4 + 0.1*sin(w1*t);
%     X_Desired(2) = 0.4 + 0.1*cos(w1*t);
%     X_Desired(3) = 0.01*sin(w2*t);
%     
%     X_Desired(4) = 0.1*w1*cos(w1*t);
%     X_Desired(5) =-0.1*w1*sin(w1*t);
%     X_Desired(6) = 0.01*w2*cos(w2*t);
%     
%     X_Desired(7) =-0.1*w1^2*sin(w1*t);
%     X_Desired(8) =-0.1*w1^2*cos(w1*t);
%     X_Desired(9) =-0.01*w2^2*sin(w2*t);

end



%% Output
Out = [a;b;c;Q_Desired;X_Desired];

end

