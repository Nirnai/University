function [X,Xcm2] = FK_Robot(u) 

    %% Forward Kinematics  
    [H0_W, H1_0,H2_0,H3_0,H4_0, Hcm2_0] = HTs(u);    
   
    H1_W = H0_W * H1_0;
    H2_W = H0_W * H2_0;
    H3_W = H0_W * H3_0;
    H4_W = H0_W * H4_0;
    
    Hcm2_W = H0_W * Hcm2_0;
    
    %% Plot Robot with current Transformation
    HT = zeros(4,4,5);
    HT(:,:,1) = H0_W;
    HT(:,:,2) = H1_W;
    HT(:,:,3) = H2_W;
    HT(:,:,4) = H3_W;
    HT(:,:,5) = H4_W;
   


    figure(1);
    clf
    view(5,15);
    robotPlot(HT);

    %% Endeffector Pose
    
    Xef = H4_W(1:3,4);
    [r,p,y] = rot2eul(H4_W(1:3,1:3));

    
    %% CM2 Pose
    Xcm2 = Hcm2_W(1:3,4);
    [r,p,y] = rot2eul(Hcm2_W(1:3,1:3));
    
    %% Return Value
    X = [Xef;Xcm2];
   
end

