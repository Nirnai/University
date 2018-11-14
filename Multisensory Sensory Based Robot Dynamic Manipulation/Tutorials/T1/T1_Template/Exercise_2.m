function Exercise_2

clc
close all

%% Provided transforms:

L=1;

l=L/(2*tan(pi/5));

H3_0=eye(4);
H3_0(1:3,4) = [0;L/2;-l];
H3_0(1:3,1:3) = RotX(-pi/10)*RotZ(-pi/2);

H4_3=eye(4);
H4_3(1:3,4) = [L*cos(pi/10);0;-L*sin(pi/10)];
H4_3(1:3,1:3) = RotY(-3*pi/10);


H3_1=eye(4);
H3_1(1:3,4) = [-(sqrt(L^2/4+(l+sqrt(l^2+L^2/4))^2))*cos(3*pi/10-atan(L/(2*(l+sqrt(l^2+(L^2)/4)))));0;-(sqrt(L^2/4+(l+sqrt(l^2+L^2/4))^2))*sin(3*pi/10-atan(L/(2*(l+sqrt(l^2+(L^2)/4)))))];
H3_1(1:3,1:3) = RotY(pi/10);

H5_2=eye(4);
H5_2(1:3,4) = [0;2*L*cos(pi/10)*sin(3*pi/10);2*L*sin(pi/10)*sin(3*pi/10)];
H5_2(1:3,1:3) = RotX(pi/5)*RotZ(3*pi/2);

H2_1=eye(4);
H2_1(1:3,4) = [-L;0;0];
H2_1(1:3,1:3) = RotZ(pi/2)*RotX(-9*pi/10);

%% Question 1: Manually computed Homogeneous Transformations

% TODO:
% Using the provided transforms, compute on matlab the absolute homogeneous
% transformations $\mathbf{H}_0^{i}$ associated with each frame of the
% pentagon. Each transform must be written as a sequence of the provided
% transforms. Put this sequence as a comment in your code.

H1_0= H3_0 * H3_1^(-1);
H2_0= H1_0 * H2_1;
H4_0= H3_0 * H4_3;
H5_0= H2_0 * H5_2;

Transformation_structure_1(:,:,1)=eye(4);
Transformation_structure_1(:,:,2)=H1_0;
Transformation_structure_1(:,:,3)=H2_0;
Transformation_structure_1(:,:,4)=H3_0;
Transformation_structure_1(:,:,5)=H4_0;
Transformation_structure_1(:,:,6)=H5_0;

figure('Name','Original Pentagon')
plotPent(Transformation_structure_1);
hold on;

%% Question 2: Relative Homogeneous Transformations

% TODO:
% Compute the relative Homogeneous Transformations between each vertex of
% the pentagon. 

H3_2 = H2_0^(-1) * H3_0;
H5_4 = H4_0^(-1) * H5_0;

%% Question 3: Relative Homogeneous Transformations

% TODO:
% Compute the transforms H5_1, H3_5 and H4_2 using at least different
% sequences of homogeneous transforms for each of them. You are free to
% chose whatever sequence you believe is right.

H5_1 = H5_0^(-1) * H1_0;
% H5_1 = H5_2^(-1) * H2_1^(-1);
% H5_1 = H5_4^(-1) * H4_3^(-1) * H3_1^(-1);

H3_5 = H5_0^(-1) * H3_0;
% H3_5 = H5_2^(-1) * H3_2; 
% H3_5 = H5_1^(-1) * H3_1;

H4_2 = H2_0^(-1) * H4_0;
% H4_2 = H5_2 * H5_4^(-1);
% H4_2 = H3_2 * H4_3;

%% Question 4: Transform the Pentagon

% TODO:

% Transform the whole pentagon including all its coordinate frames

R_a = RotX(pi/2)*RotY(pi/4)*RotZ(pi/6);
P_a = [2 0.75 0.75];

Hs_1 = eye(4);
Hs_1(1:3,1:3) = R_a;
Hs_1(1:3,4) = P_a;

Transformation_structure_2(:,:,1)=Hs_1;
Transformation_structure_2(:,:,2)=Hs_1 * H1_0;
Transformation_structure_2(:,:,3)=Hs_1 * H2_0;
Transformation_structure_2(:,:,4)=Hs_1 * H3_0;
Transformation_structure_2(:,:,5)=Hs_1 * H4_0;
Transformation_structure_2(:,:,6)=Hs_1 * H5_0;

plotPent(Transformation_structure_2);



R_b = RotZ(pi/6)*RotY(pi/4)*RotX(pi/2);
P_b = [2 -0.75 0.75];

Hs_2 = eye(4);
Hs_2(1:3,1:3) = R_b;
Hs_2(1:3,4) = P_b;

Transformation_structure_3(:,:,1)=Hs_2;
Transformation_structure_3(:,:,2)=Hs_2 * H1_0;
Transformation_structure_3(:,:,3)=Hs_2 * H2_0;
Transformation_structure_3(:,:,4)=Hs_2 * H3_0;
Transformation_structure_3(:,:,5)=Hs_2 * H4_0;
Transformation_structure_3(:,:,6)=Hs_2 * H5_0;

plotPent(Transformation_structure_3);

%% Question 5: Non-homogeneous transformations
% TODO: 

% Copy HS2
Hs_3 = Hs_2;

v = [0.5,0.5,0.5,1; 0,0,0,2; 2,2,2,1];
% Replace its last row by [a,b,c,d] where a,b,c and d are real numbers.
for j=[1,2,3]
    
    Hs_3(4,1:4) = v(j,:);
    % For each Hi_0, extract the corresponding position vector and transform it
    % using the non-homogeneous transformation.
    P1_0 = Hs_3 * H1_0(:,4);
    P1_0 = P1_0/P1_0(4);
    P2_0 = Hs_3 * H2_0(:,4);
    P2_0 = P2_0/P2_0(4);
    P3_0 = Hs_3 * H3_0(:,4);
    P3_0 = P3_0/P3_0(4);
    P4_0 = Hs_3 * H4_0(:,4);
    P4_0 = P4_0/P4_0(4);
    P5_0 = Hs_3 * H5_0(:,4);
    P5_0 = P5_0/P5_0(4);

    P = cat(2, P1_0, P2_0, P3_0, P4_0, P5_0, P1_0);
    % Plot the vertices in a new figure using the provided plot function.
    % Do this for three different vectors [a,b,c,d].

    figure();
    plot3(Hs_3(1,4), Hs_3(2,4), Hs_3(3,4),  'k .','MarkerSize',30);
    hold on;
    plot3(P(1,:),P(2,:),P(3,:),'r -', 'Linewidth',2);
end

%% Example plot function:

    function plotPent(HT)
        
        % % Visualize the cage and its respective coordinate frames
        tOff=0.1;
        aLength=0.25;
        
        grid on
        hold on
        
        axisX_0=[aLength;0;0];
        axisY_0=[0;aLength;0];
        axisZ_0=[0;0;aLength];
        
        O_0=[0;0;0];
        
        [~,~,n]=size(HT);
        
        for i=1:n
            % % Visualize the Origin Oi_0
            Oi_0(:,i)=HT(:,:,i)*[O_0;1];
            axisX1_0=HT(:,:,i)*[axisX_0;1];
            axisY1_0=HT(:,:,i)*[axisY_0;1];
            axisZ1_0=HT(:,:,i)*[axisZ_0;1];
            
            % Origin CF1
            plot3(Oi_0(1,i),Oi_0(2,i),Oi_0(3,i), 'k .','MarkerSize',30)
            text(Oi_0(1,i),Oi_0(2,i)+tOff,Oi_0(3,i)+tOff, ['O_' num2str(i-1)]);
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
        a=1.5;
        axis([-0.1 2.9 -a a -a a]);
        
        col_list='ymcrgbk';
        col=randi([1,7],1);
        random_color=col_list(col);
        final_color=['-' random_color];
        patch(Oi_0(1,2:end),Oi_0(2,2:end),Oi_0(3,2:end),final_color,'FaceAlpha',.5)
        
    end
end

