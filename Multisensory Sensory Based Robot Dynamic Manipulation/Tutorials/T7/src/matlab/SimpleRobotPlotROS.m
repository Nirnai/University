function [ Xout ] = SimpleRobotPlotROS( u )
%SIMPLEROBOTPLOT Summary of this function goes here
persistent  jointpub jointmsg counter tftree tfStampedMsg tfStampedMsg1 tfStampedMsg2 tfStampedMsg3

%Joint Position
q1=u(1);
q2=u(2);
q3=u(3);
q = [q1,q2,q3];
%Joint Velocity
qp1=u(4);
qp2=u(5);
qp3=u(6);

%Kinematic Parameters
L1=u(7);
L2=u(8);
L3=u(9);
L4=u(10);
L5=u(11);



%Time
t=u(39);

%Desired joint position
qd = [u(40),u(41),u(42)];

%Desired joint velocity (only for PID)
%qdpp = [0,0,0];

%Joint Position Vector
Q=[q1;q2;q3];

% Robot Base
T0_W = ...
[...
    1,0,0, 0;...
    0,1,0, 0;...
    0,0,1,L1;...
    0,0,0, 1 ...
];

% Homogeneous Transformations

T1_0 = ...
[...
  cos(q1), 0,  sin(q1),  0;...
  sin(q1), 0, -cos(q1),  0;...
        0, 1,        0, L1;...
        0, 0,        0,  1];

T2_0 = ...
[...
  -cos(q1)*cos(q2), cos(q1)*sin(q2),  sin(q1), -L3*cos(q1)*cos(q2);...
  -cos(q2)*sin(q1), sin(q1)*sin(q2), -cos(q1), -L3*cos(q2)*sin(q1);...
          -sin(q2),        -cos(q2),        0,     L1 - L3*sin(q2);...
                 0,               0,        0,                   1];

T3_0 = ...
[... 
  -cos(q2 + q3)*cos(q1), sin(q2 + q3)*cos(q1),  sin(q1), sin(q1)*(L2 + L4) - L3*cos(q1)*cos(q2) - L5*cos(q1)*cos(q2)*cos(q3) + L5*cos(q1)*sin(q2)*sin(q3);...
  -cos(q2 + q3)*sin(q1), sin(q2 + q3)*sin(q1), -cos(q1), L5*sin(q1)*sin(q2)*sin(q3) - L3*cos(q2)*sin(q1) - L5*cos(q2)*cos(q3)*sin(q1) - cos(q1)*(L2 + L4);...
          -sin(q2 + q3),        -cos(q2 + q3),        0,                                                                L1 - L5*sin(q2 + q3) - L3*sin(q2);...
                      0,                    0,        0,                                                                                                1];


T1_W = T0_W * T1_0;
T2_W = T0_W * T2_0;
T3_W = T0_W * T3_0;



Xef_W = T3_W(1:3,4);




if t==0
    
    %% TF publisher
    tftree = rostf;
    tfStampedMsg = rosmessage('geometry_msgs/TransformStamped');
    tfStampedMsg.Header.FrameId = 'world';
    tfStampedMsg.ChildFrameId = 'DH_0';
    
    tfStampedMsg1 = rosmessage('geometry_msgs/TransformStamped');
    tfStampedMsg1.Header.FrameId = 'ursa_base';
    tfStampedMsg1.ChildFrameId = 'DH_1';
    
    tfStampedMsg2 = rosmessage('geometry_msgs/TransformStamped');
    tfStampedMsg2.Header.FrameId = 'ursa_base';
    tfStampedMsg2.ChildFrameId = 'DH_2';
    
    tfStampedMsg3 = rosmessage('geometry_msgs/TransformStamped');
    tfStampedMsg3.Header.FrameId = 'ursa_base';
    tfStampedMsg3.ChildFrameId = 'DH_3';

    
    %% Joint State Publisher
    %Use here the correct topic name --see bringup launch file--
    jointpub = rospublisher('/ursa_joint_states', 'sensor_msgs/JointState');
    jointmsg = rosmessage(jointpub);
    
    % specific names of the joints --see urdf file--
    jointmsg.Name={ 'ursa_shoulder_pan_joint', 'ursa_shoulder_lift_joint', 'ursa_elbow_joint', 'ursa_wrist_1_joint', 'ursa_wrist_2_joint', 'ursa_wrist_3_joint'};
    
    for i=1:6
        jointmsg.Velocity(i)=0.0;
        jointmsg.Effort(i)=0.0;
    end
    
    counter=0;
    
end

%% JOINT STATE MSG and TF MSG
sampleTime=0.003;
if(~mod(t,sampleTime))
    rT=rostime('now');
    jointmsg.Header.Stamp=rT;
    jointmsg.Header.Seq=counter;
    jointmsg.Position=Q;
    send(jointpub,jointmsg);
    
    getTF(tfStampedMsg, T0_W, counter, rT);
    getTF(tfStampedMsg1, T1_0, counter, rT);
    getTF(tfStampedMsg2, T2_0, counter, rT);
    getTF(tfStampedMsg3, T3_0, counter, rT);
    
    arrayTFs=[tfStampedMsg;
        tfStampedMsg1;
        tfStampedMsg2;
        tfStampedMsg3];
    
    
    counter=counter+1;
    
    sendTransform(tftree, arrayTFs);
end


%%
Xout=Xef_W;




end

