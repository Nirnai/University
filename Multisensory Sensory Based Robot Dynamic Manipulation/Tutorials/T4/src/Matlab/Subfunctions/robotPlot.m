function [] = robotPlot(HT)

    [~,n] = size(HT);

    % Start plotting
    hold on;
    grid on;
    set(gca,'XLim',[-1 1],'YLim',[-1 1],'ZLim',[0 2])
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
    plotFrame(eye(4),'O_0');
    for i=1:n
        plotFrame(HT{2,i},['0_' num2str(i)]);
        plotFrame(HT{4,i},['0_{cm' num2str(i) '}'])
    end

 
    for i = 2:n-1
        line = [HT{2,i}(1:3,4),HT{2,i+1}(1:3,4)]';
        plot3(line(:,1), line(:,2), line(:,3),'k--');
    end


end

