function [ Xout ] = SimpleRobotPlotFull( u )

global ellipse robot

t = u(7);

%Joint Position
q1=u(4);
q2=u(5);
q3=u(6);

q = [q1;q2;q3];

%Joint Velocity
q1p=u(1);
q2p=u(2);
q3p=u(3);

qp = [q1p; q2p; q3p];


% Homogeneous Transformations H(q)
[H0_W, H1_W, H2_W, H3_W] = H_UR10_3(q);

%Compute Position of end-effector with respect to the world coordinate
Xef_0 = FK_UR10_3(q);

% Transform to wcf
Xef_W = H0_W * [Xef_0(1:3);1];

% Plot
H = zeros(4,4,4);
H(:,:,1) = H0_W;
H(:,:,2) = H1_W;
H(:,:,3) = H2_W;
H(:,:,4) = H3_W;

%Compute the Linear and Angular velocities of the end-effector wrt robot base
Xefp_0 = DifFK_UR10_3(q,qp);

%Compute the Linear velocities of the end-effector wrt wcf
Xefp_W = H0_W * [Xefp_0(1:3);1];

%% Manipulability Analisys (Ellipsoid)

% Compute Manipulability Index
Jef = Jef_0_UR10_3(u);

w = abs(det(Jef(1:3,:)));
% Compute Ellipsoid: Principal Axes 
[U,S,~] = svd(Jef(1:3,:));

A = U' * S;
vef_0= ones(4,1);

% Transform Principla axis in operational space with fixes qp = [0.1,0.1,0.1]
vef_0(1:3) = A * [0.1;0.1;0.1];


vef_W = H0_W * vef_0;

[x,y,z] = ellipsoid(Xef_W(1),Xef_W(2),Xef_W(3),vef_W(1),vef_W(2),vef_W(3));


sampleTime=0.02;

if(~mod(t,sampleTime))
    if(t == 0)
     figure(1);
    end
    clf;
    view(-60,15);
    robotPlot(H);
    plot3(Xef_W(1),Xef_W(2),Xef_W(3),'r .','MarkerSize',30 );
    text(Xef_W(1) - 0.1,Xef_W(2),Xef_W(3)+0.1, 'EF');
    ellipse = surf(x,y,z);
    alpha 0.3;
end




%% Output: vector [Xef;Xefp] (size 6X1) 
Xout=[Xef_W(1:3);Xefp_W(1:3)];

end

