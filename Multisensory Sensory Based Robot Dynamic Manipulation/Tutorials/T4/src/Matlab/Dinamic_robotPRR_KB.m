%This file was atutomatically generated by --Generate_Dinamic--
%the input vector is:
%u=[q1 q2 q3 qp4 qp5 qp6 L1 L2 L3 ]
%NOTE: The function --Genera_Robot_robotPRR_KB_Exe-- must be executed
%before running the simulink-simulator for the first time
function qpp=Dinamic_robotPRR_KB(u)

qd1 = u(4);
qd2 = u(5);
qd3 = u(6);

b1 = u(19);
b2 = u(20);
b3 = u(21);

%Joint Velocity Vector
qp=[qd1; qd2; qd3];
 
Beta = zeros(3);
Beta(1,1) = b1;
Beta(2,2) = b2;
Beta(3,3) = b3;
%Inertia Matrix +  Derivitive/Centripetal and Coriolis Matrix/Gravitational Torques Vector

[M,Md,C,G] = dynamic_model(u');
   
% Validate the skew symetric property
% N = Md -2*C;
% x = rand([3,1],'single');

% if(x'*N*x ~= 0)
%     error('C does not fullfill the skew symetrix property!');
% end
    

Tao=[0; 0; 0];

qpp=(M)\(Tao-C*qp-G-Beta*qp);

end
