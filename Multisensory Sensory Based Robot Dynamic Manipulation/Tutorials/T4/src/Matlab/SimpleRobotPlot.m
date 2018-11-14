function [ Xout ] = SimpleRobotPlot( u )
%SIMPLEROBOTPLOT Summary of this function goes here

global count
    
if(count < 4)
    count = count+1;
else   
    q1 = u(1);
    q2 = u(2);
    q3 = u(3);
    q = [q1,q2,q3];

    qd1 = u(4);
    qd2 = u(5);
    qd3 = u(6);
    qd = [qd1,qd2,qd3];


    % Robot Base

    [H1_W,H2_W,H3_W] = HTs(u');

    % 
    % H0_W = [...
    %         1  0  0 0;
    %         0 -1  0 0;
    %         0  0 -1 0;
    %         0  0  0 1];
    %    
    % H1_W = H0_W * H1_0;
    % H2_W = H0_W * H2_0;
    % H3_W = H0_W * H3_0;

    % HT = zeros(4,4,4);
    % HT(:,:,1) = H0_W;
    % HT(:,:,2) = H1_W;
    % HT(:,:,3) = H2_W;
    % HT(:,:,4) = H3_W;
    % 
    % figure(1);
    % clf
    % view(5,15);
    % robotPlot(HT);


    %Publish Joint States
    % Add your publishers here
    global tftree counter jointpub jointmsg joint1 joint2 joint3 

    if u(23) == 0
        tftree = rostf;

        joint1 = rosmessage('geometry_msgs/TransformStamped');
        joint1.ChildFrameId = 'joint_1';
        joint1.Header.FrameId = 'world';

        joint2 = rosmessage('geometry_msgs/TransformStamped');
        joint2.ChildFrameId = 'joint_2';
        joint2.Header.FrameId = 'world';

        joint3 = rosmessage('geometry_msgs/TransformStamped');
        joint3.ChildFrameId = 'joint_3';
        joint3.Header.FrameId = 'world';

        jointpub = rospublisher('/msd3dof_joint_states', 'sensor_msgs/JointState');
        jointmsg = rosmessage(jointpub);

        for i = 1:3

        jointmsg.Effort(i)=0.0;
        end

        counter = 0;
    end

    jointmsg.Position = q;
    jointmsg.Velocity = qd;

    jointmsg.Header.Stamp = rostime('now');
    jointmsg.Header.Seq=counter;
    counter=counter+1;

    showdetails(jointmsg);
    send(jointpub,jointmsg);


    joint1.Header.Stamp = rostime('now');
    joint1.Header.Seq=counter;
    joint1.Transform.Translation.X = H1_W(1,4);
    joint1.Transform.Translation.Y = H1_W(2,4);
    joint1.Transform.Translation.Z = H1_W(3,4);
    quatrot = rotm2quat(H1_W(1:3,1:3));
    joint1.Transform.Rotation.W = quatrot(1);
    joint1.Transform.Rotation.X = quatrot(2);
    joint1.Transform.Rotation.Y = quatrot(3);
    joint1.Transform.Rotation.Z = quatrot(4);


    joint2.Header.Stamp = rostime('now');
    joint2.Header.Seq=counter;
    joint2.Transform.Translation.X = H2_W(1,4);
    joint2.Transform.Translation.Y = H2_W(2,4);
    joint2.Transform.Translation.Z = H2_W(3,4);
    quatrot = rotm2quat(H2_W(1:3,1:3));
    joint2.Transform.Rotation.W = quatrot(1);
    joint2.Transform.Rotation.X = quatrot(2);
    joint2.Transform.Rotation.Y = quatrot(3);
    joint2.Transform.Rotation.Z = quatrot(4);


    joint3.Header.Stamp = rostime('now');
    joint3.Header.Seq=counter;
    joint3.Transform.Translation.X = H3_W(1,4);
    joint3.Transform.Translation.Y = H3_W(2,4);
    joint3.Transform.Translation.Z = H3_W(3,4);
    quatrot = rotm2quat(H3_W(1:3,1:3));
    joint3.Transform.Rotation.W = quatrot(1);
    joint3.Transform.Rotation.X = quatrot(2);
    joint3.Transform.Rotation.Y = quatrot(3);
    joint3.Transform.Rotation.Z = quatrot(4);


    jointpub = rospublisher('msd3dof_joint_states', 'sensor_msgs/JointState');
    jointmsg = rosmessage(jointpub);

    jointmsg.Name={ 'joint_1', ...
                    'joint_2', ...
                    'joint_3'};


    sendTransform(tftree, joint1);
    sendTransform(tftree, joint2);
    sendTransform(tftree, joint3);


    % Output: vector [Xef;Xefp;Xcm2;Xpcm2] (size 12X1) 
    
    count = 0;
end

Xout=zeros(3,1);

end

