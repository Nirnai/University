function [] = robotPlot(HT)

    [~,~,n] = size(HT);

    % Start plotting
    hold on;
    grid on;
    set(gca,'XLim',[-1 1],'YLim',[-1 1],'ZLim',[0 2])
%     ax.XLim = [-1,1];
%     ax.YLim = [-1,1];
%     ax.ZLim = [0,2];
    % Function to plot CF
        function plotFrame(H,name)
            plot3([H(1,4), H(1,4)+H(1,1)/10,H(1,4)],...
                  [H(2,4), H(2,4)+H(2,1)/10,H(2,4)],...
                  [H(3,4), H(3,4)+H(3,1)/10,H(3,4)], 'r',...
                  [H(1,4), H(1,4)+H(1,2)/10,H(1,4)],...
                  [H(2,4), H(2,4)+H(2,2)/10,H(2,4)],...
                  [H(3,4), H(3,4)+H(3,2)/10,H(3,4)], 'g',...
                  [H(1,4), H(1,4)+H(1,3)/10,H(1,4)],...
                  [H(2,4), H(2,4)+H(2,3)/10,H(2,4)],...
                  [H(3,4), H(3,4)+H(3,3)/10,H(3,4)], 'b');
             text(H(1,4),H(2,4),H(3,4),name);
        end


    % Plot CFs of 4DOF-Robot
    % plotFrame(eye(4),'O_0');


    for i=1:n
        plotFrame(HT(:,:,i),['0_' num2str(i-1)]);
    end
    
 
    for i = 2:n-1
        temp1 = HT(:,:,i);
        temp2 = HT(:,:,i+1);
        line = [temp1(1:3,4),temp2(1:3,4)]';
        plot3(line(:,1), line(:,2), line(:,3),'k--');
    end
    

end

