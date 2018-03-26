function [ Xout ] = SimpleRobotPlot( u )

%% Simulation Time

t = u(9);


%% End-effector
% Joint State
q  = [u(5);u(6);u(7);u(8)];
qp = [u(1);u(2);u(3);u(4)];


% Homogeneous Transformations wrt World frame (Numeric computation)
[H0_W, H1_W, H2_W, H3_W, H4_W] = H(q);


% Forward- and Differential Kinematics (EndEffector)
X = f(q);
Xp = J(q) * qp;

%% Centers of Mass 2
%Compute Position of the cm2 with respect to the robot base (CF 0) [numeric value]
Xcm2= fcm2(q);
Xcm2p = Jcm2(q) * qp;


%% Output

Xout=[X(1:3);Xp(1:3);Xcm2(1:3);Xcm2p(1:3)];


%% Plot

sampleTime=0.02;

if(~mod(t,sampleTime))

    HT(:,:,1) = H0_W;
    HT(:,:,2) = H1_W;
    HT(:,:,3) = H2_W;
    HT(:,:,4) = H3_W;
    HT(:,:,5) = H4_W;


    if(t == 0)    
        figure(1);
    end
    clf
    view(5,15);
    robotPlot(HT);
end

if(t == 10)
    close all;
end



end

