function Exercise_1

close all

%% Question 1: Manually computed Homogeneous Transformations

% TODO:
% By simple inspection obtain manually the absolute homogeneous transformations Ti_0 for each corner of the cage with respect to the coordinate frame O_0

T1_0=[ 1  0  0 -0.5; ...
       0  1  0 -0.4; ...
       0  0  1 -0.4; ...
       0  0  0    1];
  
T2_0=[ 0 -1  0  0.5; ...
       1  0  0 -0.4; ...
       0  0  1 -0.4; ...
       0  0  0    1];
  
T3_0=[ 0 -1  0  0.5; ...
      -1  0  0  0.4; ...
       0  0  1 -0.4; ...
       0  0  0    1];
  
T4_0=[ 0  1  0 -0.5; ...
      -1  0  0  0.4; ...
       0  0  1 -0.4; ...
       0  0  0    1];
  
T5_0=[ 1  0  0 -0.5; ...
       0  1  0 -0.4; ...
       0  0 -1  0.4; ...
       0  0  0    1];
  
T6_0=[-1  0  0  0.5; ...
       0  1  0 -0.4; ...
       0  0 -1  0.4; ...
       0  0  0    1];
  
T7_0=[ 0 -1  0  0.5; ...
      -1  0  0  0.4; ...
       0  0 -1  0.4; ...
       0  0  0    1];
  
T8_0=[ 1  0  0 -0.5; ...
       0 -1  0  0.4; ...
       0  0 -1  0.4; ...
       0  0  0    1];

Transformation_structure_1(:,:,1)=T1_0;
Transformation_structure_1(:,:,2)=T2_0;
Transformation_structure_1(:,:,3)=T3_0;
Transformation_structure_1(:,:,4)=T4_0;
Transformation_structure_1(:,:,5)=T5_0;
Transformation_structure_1(:,:,6)=T6_0;
Transformation_structure_1(:,:,7)=T7_0;
Transformation_structure_1(:,:,8)=T8_0;

%% Question 2: Left handed frames

% TODO:
%  The rotation matrix of some frames in the previous figure can not be 
% obtained using a sequence of basic rotations. Please explain why 
% (elaborate your answer using for example the definition of SE(3)), and 
% propose the proper coordinate frames to fix this problem.


figure('Name','Cage with Left handed O3 and O5 (in Magenta)')
plotCage(Transformation_structure_1);
hold on
plot3(0.5,0.4,-0.4,'m .','MarkerSize',60);
plot3(-0.5,-0.4, 0.4,'m .','MarkerSize',60);

% Explanation: The initial coordinate frame is a right handed frame. It is
% not possible to transform a right handed frame into a left handed on. The
% rotation matrix necessary for that has a determinant of -1 und thus is
% not part of SE(3). 

%% Question 3: Right Hand Frames

% TODO:
%  Replace the coordinate frames found in the above question with the 
% suitable coordinate frames, respecting the figure convention; 
% then compute the associated absolute homogeneous transformations Hi_0

H1_0=[ 1  0  0 -0.5; ...
       0  1  0 -0.4; ...
       0  0  1 -0.4; ...
       0  0  0    1];
  
H2_0=[ 0 -1  0  0.5; ...
       1  0  0 -0.4; ...
       0  0  1 -0.4; ...
       0  0  0    1];
  
H3_0=[ 0  1  0  0.5; ...
      -1  0  0  0.4; ...
       0  0  1 -0.4; ...
       0  0  0    1];
  
H4_0=[ 0  1  0 -0.5; ...
      -1  0  0  0.4; ...
       0  0  1 -0.4; ...
       0  0  0    1];
  
H5_0=[ 1  0  0 -0.5; ...
       0 -1  0 -0.4; ...
       0  0 -1  0.4; ...
       0  0  0    1];
  
H6_0=[-1  0  0  0.5; ...
       0  1  0 -0.4; ...
       0  0 -1  0.4; ...
       0  0  0    1];
  
H7_0=[ 0 -1  0  0.5; ...
      -1  0  0  0.4; ...
       0  0 -1  0.4; ...
       0  0  0    1];
  
H8_0=[ 1  0  0 -0.5; ...
       0 -1  0  0.4; ...
       0  0 -1  0.4; ...
       0  0  0    1];

Transformation_structure_2(:,:,1)=H1_0;
Transformation_structure_2(:,:,2)=H2_0;
Transformation_structure_2(:,:,3)=H3_0;
Transformation_structure_2(:,:,4)=H4_0;
Transformation_structure_2(:,:,5)=H5_0;
Transformation_structure_2(:,:,6)=H6_0;
Transformation_structure_2(:,:,7)=H7_0;
Transformation_structure_2(:,:,8)=H8_0;

figure('Name','Correct Cage')
plotCage(Transformation_structure_2);

%% Question 4: Rotation matrices and translation vectors

% TODO:
% Extract the rotation matrices and translation vectors from each HT:
% 
R1_0 = [ 1  0  0; ...
         0  1  0; ...
         0  0  1];
   
P1_0 = [-0.5 -0.4 -0.4];
  


R2_0 = [ 0 -1  0; ...
         1  0  0; ...
         0  0  1];
    
P2_0 = [0.5 -0.4 -0.4];
  


R3_0 = [ 0  1  0; ...
        -1  0  0; ...
         0  0  1];
     
P3_0 = [0.5 0.4 -0.4];
  


R4_0 = [ 0  1  0; ...
        -1  0  0; ...
         0  0  1];
   
P4_0 = [-0.5 0.4 -0.4];
  


R5_0 = [ 1  0  0; ...
         0 -1  0; ...
         0  0 -1];
   
P5_0 = [-0.5 -0.4 0.4];


  
R6_0 = [-1  0  0; ...
         0  1  0; ...
         0  0 -1];

P6_0 = [0.5 -0.4 0.4];

  


R7_0 = [ 0 -1  0; ...
        -1  0  0; ...
         0  0 -1];
   
P7_0 = [0.5 0.4 0.4];


  
R8_0 = [ 1  0  0; ...
         0 -1  0; ...
         0  0 -1];
   
P8_0 = [ -0.5 0.4 0.4];    


%% Question 5: Simple Euler angles without considering singularities

% TODO:
% Represent each rotation matrix as a set of Euler angles (roll-pitch-yaw)

[phi_1,theta_1,psi_1]=R2EulerA(R1_0);
[phi_2,theta_2,psi_2]=R2EulerA(R2_0);
[phi_3,theta_3,psi_3]=R2EulerA(R3_0);
[phi_4,theta_4,psi_4]=R2EulerA(R4_0);
[phi_5,theta_5,psi_5]=R2EulerA(R5_0);
[phi_6,theta_6,psi_6]=R2EulerA(R6_0);
[phi_7,theta_7,psi_7]=R2EulerA(R7_0);
[phi_8,theta_8,psi_8]=R2EulerA(R8_0);



%% Question 6: Relative Homogeneous Transformations (manually Generated)

% TODO:
% Compute manually the the relative transformations between the corners of
% the cage

H0_1 = H1_0^(-1);
H2_1= H0_1 * H2_0;

H0_2 = H2_0^(-1);
H3_2= H0_2 * H3_0;

H0_3 = H3_0^(-1);
H4_3= H0_3 * H4_0;

H0_4 = H4_0^(-1);
H5_4= H0_4 * H5_0;

H0_5 = H5_0^(-1);
H6_5= H0_5 * H6_0;

H0_6 = H6_0^(-1);
H7_6= H0_6 * H7_0;

H0_7 = H7_0^(-1);
H8_7= H0_7 * H8_0;

%% Question 7: Absolute Transfomations using the relative transformations

% TODO:
% Compute again the absolute Homogeneous Transformations of each coordinate
% frame using the relative transformations and compare them with the

H2_0 = H1_0 * H2_1;
H3_0 = H1_0 * H2_1 * H3_2;
H4_0 = H1_0 * H2_1 * H3_2 * H4_3;
H5_0 = H1_0 * H2_1 * H3_2 * H4_3 * H5_4;
H6_0 = H1_0 * H2_1 * H3_2 * H4_3 * H5_4 * H6_5;
H7_0 = H1_0 * H2_1 * H3_2 * H4_3 * H5_4 * H6_5 * H7_6;
H8_0 = H1_0 * H2_1 * H3_2 * H4_3 * H5_4 * H6_5 * H7_6 * H8_7;

Transformation_structure_3(:,:,1)=H1_0;
Transformation_structure_3(:,:,2)=H2_0;
Transformation_structure_3(:,:,3)=H3_0;
Transformation_structure_3(:,:,4)=H4_0;
Transformation_structure_3(:,:,5)=H5_0;
Transformation_structure_3(:,:,6)=H6_0;
Transformation_structure_3(:,:,7)=H7_0;
Transformation_structure_3(:,:,8)=H8_0;

figure('Name','Correct Cage with relative HT')
plotCage(Transformation_structure_3);

%% Question 8: Plot the Cage

% TODO: Plot each coordinate frame in Matlab to visualize the cage

    function plotCage(HT)
        
        % % Visualize the cage and its respective coordinate frames
        tOff=0.05;
        aLength=0.25;
        
        grid on
        hold on
        
        axisX_0=[aLength;0;0];
        axisY_0=[0;aLength;0];
        axisZ_0=[0;0;aLength];
        
        % % Visualize the Origin O_0
        
        O_0=[0;0;0];
        
        plot3(O_0(1),O_0(2),O_0(3), 'k .','MarkerSize',30)
        text(O_0(1)+tOff,O_0(2)+tOff,O_0(3)+tOff, 'O_0');
        %Plot x-axis
        plot3([O_0(1);axisX_0(1)],[O_0(2);axisX_0(2)],[O_0(3);axisX_0(3)],'r -', 'Linewidth',2)
        %Plot y-axis
        plot3([O_0(1);axisY_0(1)],[O_0(2);axisY_0(2)],[O_0(3);axisY_0(3)],'g -', 'Linewidth',2)
        %Plot z-axis
        plot3([O_0(1);axisZ_0(1)],[O_0(2);axisZ_0(2)],[O_0(3);axisZ_0(3)],'b -', 'Linewidth',2)
        
        [~,~,n]=size(HT);
        
        for i=1:n
            % % Visualize the Origin Oi_0
            Oi_0(:,i)=HT(:,:,i)*[O_0;1];
            axisX1_0=HT(:,:,i)*[axisX_0;1];
            axisY1_0=HT(:,:,i)*[axisY_0;1];
            axisZ1_0=HT(:,:,i)*[axisZ_0;1];
            
            % Origin CF1
            plot3(Oi_0(1,i),Oi_0(2,i),Oi_0(3,i), 'k .','MarkerSize',30)
            text(Oi_0(1,i)+tOff,Oi_0(2,i)+tOff,Oi_0(3,i)+tOff, ['O_' num2str(i)]);
            %Plot x-axis
            plot3([Oi_0(1,i);axisX1_0(1)],[Oi_0(2,i);axisX1_0(2)],[Oi_0(3,i);axisX1_0(3)],'r -', 'Linewidth',2)
            %Plot y-axis
            plot3([Oi_0(1,i);axisY1_0(1)],[Oi_0(2,i);axisY1_0(2)],[Oi_0(3,i);axisY1_0(3)],'g -', 'Linewidth',2)
            %Plot z-axis
            plot3([Oi_0(1,i);axisZ1_0(1)],[Oi_0(2,i);axisZ1_0(2)],[Oi_0(3,i);axisZ1_0(3)],'b -', 'Linewidth',2)
        end
        
        xlabel('X');
        ylabel('Y');
        zlabel('Z');
        
        % Cage structure:
        
        plot3([Oi_0(1,1);Oi_0(1,2)],[Oi_0(2,1);Oi_0(2,2)],[Oi_0(3,1);Oi_0(3,2)],'k --', 'Linewidth',0.5)
        plot3([Oi_0(1,1);Oi_0(1,5)],[Oi_0(2,1);Oi_0(2,5)],[Oi_0(3,1);Oi_0(3,5)],'k --', 'Linewidth',0.5)
        plot3([Oi_0(1,2);Oi_0(1,3)],[Oi_0(2,2);Oi_0(2,3)],[Oi_0(3,2);Oi_0(3,3)],'k --', 'Linewidth',0.5)
        plot3([Oi_0(1,2);Oi_0(1,6)],[Oi_0(2,2);Oi_0(2,6)],[Oi_0(3,2);Oi_0(3,6)],'k --', 'Linewidth',0.5)
        plot3([Oi_0(1,3);Oi_0(1,4)],[Oi_0(2,3);Oi_0(2,4)],[Oi_0(3,3);Oi_0(3,4)],'k --', 'Linewidth',0.5)
        plot3([Oi_0(1,3);Oi_0(1,7)],[Oi_0(2,3);Oi_0(2,7)],[Oi_0(3,3);Oi_0(3,7)],'k --', 'Linewidth',0.5)
        plot3([Oi_0(1,4);Oi_0(1,8)],[Oi_0(2,4);Oi_0(2,8)],[Oi_0(3,4);Oi_0(3,8)],'k --', 'Linewidth',0.5)
        plot3([Oi_0(1,4);Oi_0(1,1)],[Oi_0(2,4);Oi_0(2,1)],[Oi_0(3,4);Oi_0(3,1)],'k --', 'Linewidth',0.5)
        plot3([Oi_0(1,5);Oi_0(1,6)],[Oi_0(2,5);Oi_0(2,6)],[Oi_0(3,5);Oi_0(3,6)],'k --', 'Linewidth',0.5)
        plot3([Oi_0(1,6);Oi_0(1,7)],[Oi_0(2,6);Oi_0(2,7)],[Oi_0(3,6);Oi_0(3,7)],'k --', 'Linewidth',0.5)
        plot3([Oi_0(1,7);Oi_0(1,8)],[Oi_0(2,7);Oi_0(2,8)],[Oi_0(3,7);Oi_0(3,8)],'k --', 'Linewidth',0.5)
        plot3([Oi_0(1,8);Oi_0(1,5)],[Oi_0(2,8);Oi_0(2,5)],[Oi_0(3,8);Oi_0(3,5)],'k --', 'Linewidth',0.5)
        
    end
end

