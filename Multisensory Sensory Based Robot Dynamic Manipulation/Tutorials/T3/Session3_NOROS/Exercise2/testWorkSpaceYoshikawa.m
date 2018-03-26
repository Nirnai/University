function testWorkSpaceYoshikawa( fig,samples,wThreshold, sORw, plotRobot)
%TESTWORKSPACE Summary of this function goes here
% Workspace and Singularities analysis
% Arguments:
% fig: Id for the figure
% samples: Number of max samples to be tested
% wThreshold: threshold for the determinant.
% sORw: Singularity or Workspace analysis true= Sin, false= WS
% PlotRobot: Flag to define if the robot will be shown in the figure.


close all

% Define the figure and its properties
figure(fig)
hold on
workSpace=[-1.5 1.5 -1.5 1.5 -1.5 1.5];
axis(workSpace);

% Color map to plot the different manipulability indexes 
nColors=1000;
cmap=colormap(jet(nColors));


q1=zeros(samples,1);
q2=zeros(samples,1);
q3=zeros(samples,1);

w=zeros(samples,1);

q1min=deg2rad(-180);
q1max=deg2rad(180);

q2min=deg2rad(-180);
q2max=deg2rad(180);

q3min=deg2rad(-180);
q3max=deg2rad(180);



% Compute Manipulability indexes for each sample data
parfor i=1:samples
    % Random joint pose     
    q1(i)=((q1max-q1min)*rand)+q1min;
    q2(i)=((q2max-q2min)*rand)+q2min;
    q3(i)=((q3max-q3min)*rand)+q3min;
     
    Jef = Jef_0_UR10_3([q1(i);q2(i);q3(i)]);

    w(i)=det(Jef*Jef');
end


% Get the max index value
wMax=max(w);

xef_0=zeros(6,samples);
colorIDX_0=[];

% Discriminate all the Samples for Singularity and WS
if(sORw)
    parfor i=1:samples
        
        if w(i)<wThreshold
            % store the end-effector position in a Matrix
            xef_0(:,i) = FK_UR10_3([q1(i);q2(i);q3(i)]);
            % Store the corresponding color in an array (it depends on wMax)
            % if w(i)==wMax -> Dark Blue; if w(i)== Dark Red
           if w(i)==wMax 
            colorIDX_0(i) ="DarkBlue";
           else
            colorIDX_0(i) ="DarkRed";
           end
            
        end
    end
else
    parfor i=1:samples
        if w(i)>wThreshold
	    % store the end-effector position in a Matrix
             xef_0(:,i) = FK_UR10_3([q1(i);q2(i);q3(i)]);
            % Store the corresponding color in an array (it depends on wMax)
	    % if w(i)==wMax -> Dark Blue; if w(i)== Dark Red	
           if (w(i)==wMax) 
            colorIDX_0(i) ="DarkBlue";
           else
            colorIDX_0(i) ="DarkRed";
           end
        end
    end
end

% Number of computed points
sxef=size(xef_0);

% Plot Robot for visual reference
if(plotRobot)
    % Plot the robot in a Joint position with the max manipulability 
    
end

% Robot Base
T0_W=eye(4);


% Plot points
for j=1:sxef(2)
    plot3(xef_0(1,j),xef_0(2,j),xef_0(3,j),'*','Color',cmap(colorIDX_0(j,1),:)), view(-60,15);
end

end

