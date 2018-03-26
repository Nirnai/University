function [ Xout ] = SimpleRobotPlotROS( u )

global tftree counter jointpub jointmsg joint1 joint2 joint3 EF
    
  
    q1 = u(1);
    q2 = u(2);
    q3 = u(3);
    q = [q1,q2,q3];
    
    H0_W = [...
        1,0,0,     0;...
        0,1,0,     0;...
        0,0,1,0.1273;...
        0,0,0,     1];

    [H1_0,H2_0,H3_0] = H(q);


    %Initialize the publishers and messages
    
    if u(7) == 0
        tftree = rostf;

        joint1 = rosmessage('geometry_msgs/TransformStamped');
        joint1.ChildFrameId = 'ursa_shoulder_pan_joint';
        joint1.Header.FrameId = 'world';

        joint2 = rosmessage('geometry_msgs/TransformStamped');
        joint2.ChildFrameId = 'ursa_shoulder_lift_joint';
        joint2.Header.FrameId = 'world';

        joint3 = rosmessage('geometry_msgs/TransformStamped');
        joint3.ChildFrameId = 'ursa_elbow_joint';
        joint3.Header.FrameId = 'world';
        
        EF = rosmessage('geometry_msgs/TransformStamped');
        EF.ChildFrameId = 'EF';
        EF.Header.FrameId = 'world';

        jointpub = rospublisher('/ursa_joint_states', 'sensor_msgs/JointState');
        jointmsg = rosmessage(jointpub);

        counter = 0;
    end
    
    
    sampleTime=0.02;
    if(~mod(u(7),sampleTime))
        for i = 1:3
        jointmsg.Effort(i)=0.0;
        jointmsg.Velocity(i) = 0.0;
        end

        jointmsg.Position = q;


        jointmsg.Header.Stamp = rostime('now');
        jointmsg.Header.Seq=counter;
        counter=counter+1;

        send(jointpub,jointmsg);

        %Publish the robot joints and tf's (Base, Links, CMs and EF), see the joint names generated by bringUR10.launch
        joint1.Header.Stamp = rostime('now');
        joint1.Header.Seq=counter;
        joint1.Transform.Translation.X = H0_W(1,4);
        joint1.Transform.Translation.Y = H0_W(2,4);
        joint1.Transform.Translation.Z = H0_W(3,4);
        quatrot = rotm2quat(H0_W(1:3,1:3));
        joint1.Transform.Rotation.W = quatrot(1);
        joint1.Transform.Rotation.X = quatrot(2);
        joint1.Transform.Rotation.Y = quatrot(3);
        joint1.Transform.Rotation.Z = quatrot(4);

        joint2.Header.Stamp = rostime('now');
        joint2.Header.Seq=counter;
        joint2.Transform.Translation.X = H1_0(1,4);
        joint2.Transform.Translation.Y = H1_0(2,4);
        joint2.Transform.Translation.Z = H1_0(3,4);
        quatrot = rotm2quat(H1_0(1:3,1:3));
        joint2.Transform.Rotation.W = quatrot(1);
        joint2.Transform.Rotation.X = quatrot(2);
        joint2.Transform.Rotation.Y = quatrot(3);
        joint2.Transform.Rotation.Z = quatrot(4);


        joint3.Header.Stamp = rostime('now');
        joint3.Header.Seq=counter;
        joint3.Transform.Translation.X = H2_0(1,4);
        joint3.Transform.Translation.Y = H2_0(2,4);
        joint3.Transform.Translation.Z = H2_0(3,4);
        quatrot = rotm2quat(H2_0(1:3,1:3));
        joint3.Transform.Rotation.W = quatrot(1);
        joint3.Transform.Rotation.X = quatrot(2);
        joint3.Transform.Rotation.Y = quatrot(3);
        joint3.Transform.Rotation.Z = quatrot(4);

        EF.Header.Stamp = rostime('now');
        EF.Header.Seq=counter;
        EF.Transform.Translation.X = H3_0(1,4);
        EF.Transform.Translation.Y = H3_0(2,4);
        EF.Transform.Translation.Z = H3_0(3,4);
        quatrot = rotm2quat(H3_0(1:3,1:3));
        EF.Transform.Rotation.W = quatrot(1);
        EF.Transform.Rotation.X = quatrot(2);
        EF.Transform.Rotation.Y = quatrot(3);
        EF.Transform.Rotation.Z = quatrot(4);


        jointpub = rospublisher('ursa_joint_states', 'sensor_msgs/JointState');
        jointmsg = rosmessage(jointpub);

        jointmsg.Name={ 'ursa_shoulder_pan_joint', ...
                        'ursa_shoulder_lift_joint', ...
                        'ursa_elbow_joint'};


        sendTransform(tftree, joint1);
        sendTransform(tftree, joint2);
        sendTransform(tftree, joint3);
        sendTransform(tftree, EF);
    end

% Output: vector [Xef] (size 3X1)
Xout=[0;0;0];

end
