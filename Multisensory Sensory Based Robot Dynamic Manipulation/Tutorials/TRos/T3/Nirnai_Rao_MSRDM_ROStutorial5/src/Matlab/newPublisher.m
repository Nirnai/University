function x = newPublisher( u )

    global  jointpub jointmsg counter tftree q1Joint q2Joint q3Joint q4Joint cm1 cm2 cm3 cm4
    if u==0 
        
        %% TF publisher
        tftree = rostf;
        
        q1Joint = rosmessage('geometry_msgs/TransformStamped');
        q1Joint.ChildFrameId = 'q1_joint';
        q1Joint.Header.FrameId = 'world';
        
        q2Joint = rosmessage('geometry_msgs/TransformStamped');
        q2Joint.ChildFrameId = 'q2_joint';
        q2Joint.Header.FrameId = 'world';
        
        q3Joint = rosmessage('geometry_msgs/TransformStamped');
        q3Joint.ChildFrameId = 'q3_joint';
        q3Joint.Header.FrameId = 'world';
        
        q4Joint = rosmessage('geometry_msgs/TransformStamped');
        q4Joint.ChildFrameId = 'q4_joint';
        q4Joint.Header.FrameId = 'world';
        
        cm1 = rosmessage('geometry_msgs/TransformStamped');
        cm1.ChildFrameId = 'cm1';
        cm1.Header.FrameId = 'world';
        
        cm1 = rosmessage('geometry_msgs/TransformStamped');
        cm1.ChildFrameId = 'cm1';
        cm1.Header.FrameId = 'world';
        
        cm2 = rosmessage('geometry_msgs/TransformStamped');
        cm2.ChildFrameId = 'cm2';
        cm2.Header.FrameId = 'world';
        
        cm3 = rosmessage('geometry_msgs/TransformStamped');
        cm3.ChildFrameId = 'cm3';
        cm3.Header.FrameId = 'world';
        
        cm4 = rosmessage('geometry_msgs/TransformStamped');
        cm4.ChildFrameId = 'cm4';
        cm4.Header.FrameId = 'world';
        %% Joint State Publisher
        %Use here the correct topic name --see bringup launch file--
        jointpub = rospublisher('/fourdof_joint_states', 'sensor_msgs/JointState');
        jointmsg = rosmessage(jointpub);

        % specific names of the joints --see urdf file--
        jointmsg.Name={ 'q1_joint', ...
                        'q2_joint', ...
                        'q3_joint', ...
                        'q4_joint', ...
                        'cm1', ...
                        'cm2', ...
                        'cm3', ...
                        'cm4'};

        for i=1:8
        jointmsg.Velocity(i)=0.0;
        jointmsg.Effort(i)=0.0;
        end

        counter=0;

    end
    
    %% JOINT STATE MSG
    T1 = 2;
    T2 = 2;
    T3 = 2;
    T4 = 2;
    w1 = 2*pi/T1;
    v2 = 2*pi/T2;
    v3 = 2*pi/T3;
    w4 = 2*pi/T4;


    jointmsg.Header.Stamp = rostime('now');
    jointmsg.Header.Seq=counter;
    counter=counter+1;

    %Use the joint limits of the specific robot --see urdf file--
    q1 = pi/2*sin(w1*u);
    q2 = 0.06 + 0.04*sin(v2*u);
    q3 = 0.06 + 0.04*sin(v3*u);
    q4 = pi/2*sin(w4*u);
    jointmsg.Position=[q1,q2,q3,q4];
    showdetails(jointmsg);
    send(jointpub,jointmsg);

    %% D_H Transform

    % Kinematic Parameters
    % add as many variable as needed, according to the robot parameters
    L1=0.30;
    L2=0.34;
    L3=0.583;
    L4=0.08+0.02;
    L5=0.08;
    L6=0.115;
    L7=0.082;
    L8=0.108;
    L9=0.123;
    L10=0.167;
    L11=0.028;
    L12=0.050;
    c = 0.32;           % Table length
    T1 = 45*pi/180;     % Theta1
    T2 = 30*pi/180;     % Theat2

    % compute some intermediate lengths
    L13 = L9*sin(T2);
    L14 = L9*cos(T2);
    L15 = (L12 + c/2)*cos(T1);
    L16 = (L12 + c/2)*sin(T1);
    L17 = sqrt((cos(T2)*L9)^2 + (L10/cos(T1))^2);

    % compute intermediate angle
    Anlge = atan(cos(T2)*cos(T1)*L9/L10);

    % combine length to simplify symbolic equations
    L18 = L6 + q2 + 0.1 + L5;
    L19 = L7 + q3 + 0.1 + L8 + L13;
    L20 = L10 + L11 + L15;
    L21 = L6/2 + q2 + 0.1 + L5;

    %% Joints and End-effector
    % Robot Base SEE FIG.6 
    T0_W=[0 -1 0 L1; 0 0 -1 L2; 1 0 0 L3; 0 0 0 1];

    % Relative Homogeneous Transformations (symbolic form)

    T1_0 = [[cos(q1), 0, -sin(q1), 0]; [sin(q1), 0, cos(q1), 0]; [0, -1, 0, L4]; [0, 0, 0, 1]];

    T2_1 = [[0, 0, -1, 0]; [-1, 0, 0, 0]; [0, 1, 0, L18]; [0, 0, 0, 1]];

    T3_2 = [[1, 0, 0, L14]; [0, -2^(1/2)/2, -2^(1/2)/2, 0]; [0, 2^(1/2)/2, -2^(1/2)/2, L19]; [0, 0, 0, 1]];

    T4_3 = [[-sin(q4), -cos(q4), 0, -L16*sin(q4)]; [cos(q4), -sin(q4), 0, L16*cos(q4)]; [0, 0, 1, L20]; [0, 0, 0, 1]];

    % Homogeneous Transformations wrt BASE frame (Numeric computation)

    T2_0 = T1_0 * T2_1;
    T3_0 = T1_0 * T2_1 * T3_2;
    T4_0 = T1_0 * T2_1 * T3_2 * T4_3;

    % Homogeneous Transformations wrt WORLD frame (Numeric computation)

    T1_W = T0_W * T1_0;
    T2_W = T0_W * T2_0;
    T3_W = T0_W * T3_0;
    T4_W = T0_W * T4_0;
    
    t1_W = T1_W(1:3,4);
    t2_W = T2_W(1:3,4);
    t3_W = T3_W(1:3,4);
    t4_W = T4_W(1:3,4);
    
    R1_W = T1_W(1:3,1:3);
    R2_W = T2_W(1:3,1:3);
    R3_W = T3_W(1:3,1:3);
    R4_W = T4_W(1:3,1:3);
    
    % Center of Masses
    Tcm1_0 = [[-sin(q1), -cos(q1), 0, -(L5*sin(q1))/2]; [cos(q1), -sin(q1), 0, (L5*cos(q1))/2]; [0, 0, 1, L4]; [0, 0, 0, 1]];
    Tcm2_1 = [[1, 0, 0, -L7/2]; [0, 1, 0, 0]; [0, 0, 1, L21]; [0, 0, 0, 1]];
    Tcm3_2 = [[-sin(Anlge), -cos(Anlge), 0, L17*sin(Anlge)]; [cos(Anlge), -sin(Anlge), 0, -L17*cos(Anlge)]; [0, 0, 1, L19]; [0, 0, 0, 1]];
    Tcm4_3 = [[-sin(q4), -cos(q4), 0, -L16*sin(q4)]; [cos(q4), -sin(q4), 0, L16*cos(q4)]; [0, 0, 1, L20]; [0, 0, 0, 1]];

    % Homogeneous Transformations wrt base frame (Numeric computation)
    Tcm2_0 = T1_0 * Tcm2_1;
    Tcm3_0 = T2_0 * Tcm3_2;
    Tcm4_0 = T3_0 * Tcm4_3;

    % Homogeneous Transformations wrt World frame (Numeric computation)
    Tcm1_W = T0_W * Tcm1_0;
    Tcm2_W = T0_W * Tcm2_0;
    Tcm3_W = T0_W * Tcm3_0;
    Tcm4_W = T0_W * Tcm4_0;
    
    tcm1_W = Tcm1_W(1:3,4);
    tcm2_W = Tcm2_W(1:3,4);
    tcm3_W = Tcm3_W(1:3,4);
    tcm4_W = Tcm4_W(1:3,4);
    
    Rcm1_W = Tcm1_W(1:3,1:3);
    Rcm2_W = Tcm2_W(1:3,1:3);
    Rcm3_W = Tcm3_W(1:3,1:3);
    Rcm4_W = Tcm4_W(1:3,1:3);
    
    %% TF MSG

    q1Joint.Header.Stamp = rostime('now');
    q1Joint.Header.Seq=counter;
    q1Joint.Transform.Translation.X = t1_W(1);
    q1Joint.Transform.Translation.Y = t1_W(2);
    q1Joint.Transform.Translation.Z = t1_W(3);
    quatrot = rotm2quat(R1_W);
    q1Joint.Transform.Rotation.W = quatrot(1);
    q1Joint.Transform.Rotation.X = quatrot(2);
    q1Joint.Transform.Rotation.Y = quatrot(3);
    q1Joint.Transform.Rotation.Z = quatrot(4);

    q2Joint.Header.Stamp = rostime('now');
    q2Joint.Header.Seq=counter;
    q2Joint.Transform.Translation.X = t2_W(1);
    q2Joint.Transform.Translation.Y = t2_W(2);
    q2Joint.Transform.Translation.Z = t2_W(3);
    quatrot = rotm2quat(R2_W);
    q2Joint.Transform.Rotation.W = quatrot(1);
    q2Joint.Transform.Rotation.X = quatrot(2);
    q2Joint.Transform.Rotation.Y = quatrot(3);
    q2Joint.Transform.Rotation.Z = quatrot(4);

    q3Joint.Header.Stamp = rostime('now');
    q3Joint.Header.Seq=counter;
    q3Joint.Transform.Translation.X = t3_W(1);
    q3Joint.Transform.Translation.Y = t3_W(2);
    q3Joint.Transform.Translation.Z = t3_W(3);
    quatrot = rotm2quat(R3_W);
    q3Joint.Transform.Rotation.W = quatrot(1);
    q3Joint.Transform.Rotation.X = quatrot(2);
    q3Joint.Transform.Rotation.Y = quatrot(3);
    q3Joint.Transform.Rotation.Z = quatrot(4);

    q4Joint.Header.Stamp = rostime('now');
    q4Joint.Header.Seq=counter;
    q4Joint.Transform.Translation.X = t4_W(1);
    q4Joint.Transform.Translation.Y = t4_W(2);
    q4Joint.Transform.Translation.Z = t4_W(3);
    quatrot = rotm2quat(R4_W);
    q4Joint.Transform.Rotation.W = quatrot(1);
    q4Joint.Transform.Rotation.X = quatrot(2);
    q4Joint.Transform.Rotation.Y = quatrot(3);
    q4Joint.Transform.Rotation.Z = quatrot(4);

    cm1.Header.Stamp = rostime('now');
    cm1.Header.Seq=counter;
    cm1.Transform.Translation.X = tcm1_W(1);
    cm1.Transform.Translation.Y = tcm1_W(2);
    cm1.Transform.Translation.Z = tcm1_W(3);
    quatrot = rotm2quat(Rcm1_W);
    cm1.Transform.Rotation.W = quatrot(1);
    cm1.Transform.Rotation.X = quatrot(2);
    cm1.Transform.Rotation.Y = quatrot(3);
    cm1.Transform.Rotation.Z = quatrot(4);

    cm2.Header.Stamp = rostime('now');
    cm2.Header.Seq=counter;
    cm2.Transform.Translation.X = tcm2_W(1);
    cm2.Transform.Translation.Y = tcm2_W(2);
    cm2.Transform.Translation.Z = tcm2_W(3);
    quatrot = rotm2quat(Rcm2_W);
    cm2.Transform.Rotation.W = quatrot(1);
    cm2.Transform.Rotation.X = quatrot(2);
    cm2.Transform.Rotation.Y = quatrot(3);
    cm2.Transform.Rotation.Z = quatrot(4);

    cm3.Header.Stamp = rostime('now');
    cm3.Header.Seq=counter;
    cm3.Transform.Translation.X = tcm3_W(1);
    cm3.Transform.Translation.Y = tcm3_W(2);
    cm3.Transform.Translation.Z = tcm3_W(3);
    quatrot = rotm2quat(Rcm3_W);
    cm3.Transform.Rotation.W = quatrot(1);
    cm3.Transform.Rotation.X = quatrot(2);
    cm3.Transform.Rotation.Y = quatrot(3);
    cm3.Transform.Rotation.Z = quatrot(4);

    cm4.Header.Stamp = rostime('now');
    cm4.Header.Seq=counter;
    cm4.Transform.Translation.X = tcm4_W(1);
    cm4.Transform.Translation.Y = tcm4_W(2);
    cm4.Transform.Translation.Z = tcm4_W(3);
    quatrot = rotm2quat(Rcm4_W);
    cm4.Transform.Rotation.W = quatrot(1);
    cm4.Transform.Rotation.X = quatrot(2);
    cm4.Transform.Rotation.Y = quatrot(3);
    cm4.Transform.Rotation.Z = quatrot(4);



    sendTransform(tftree, q1Joint)
    sendTransform(tftree, q2Joint)
    sendTransform(tftree, q3Joint)
    sendTransform(tftree, q4Joint)
    sendTransform(tftree, cm1)
    sendTransform(tftree, cm2)
    sendTransform(tftree, cm3)
    sendTransform(tftree, cm4)

    %Use the above example to plot the TF of each joint and each center 
    %of mass using the Homogenous transformations obtained from 
    %D-H parameters. You will need the joint positions which are generated
    %with the JOINT STATE MSG publisher (see above)




    x=u;

end

