% Function:
% Visualize the distribution of sensor nodes before and after optimization
%
% Inputs:
% Positions_Before, Positions_After: Sensor positions before and after optimization
% day: Current day of optimization
% longitude_range, latitude_range, depth_range: Ranges of longitude, latitude, and depth
% centerPosition_noSOA, centerPosition_SOA: Positions of central nodes before and after optimization
%
% Outputs:
% None
%
function VisualizeSensorPositions(Positions_Before, Positions_After, day, longitude_range, latitude_range, depth_range, centerPosition_noSOA, centerPosition_SOA)
    % Plot the 3D scatter plot of sensor positions before optimization
    fig1 = figure;
    hold on;
    scatter3(Positions_Before(:,1), Positions_Before(:,2), Positions_Before(:,3), 'blue', 'filled');
    scatter3(centerPosition_noSOA(:,1), centerPosition_noSOA(:,2), centerPosition_noSOA(:,3), 'red', 'filled');
    titleText = ['Day ', num2str(day), ' Sensor Positions Before Optimization'];
    title(titleText);
    xlabel('Longitude(m)');
    ylabel('Latitude(m)');
    zlabel('Depth(m)');
    set(gca, 'ZDir', 'reverse'); % Reverse the Z axis so depth increases downwards
    axis([longitude_range(1) longitude_range(2) latitude_range(1) latitude_range(2) depth_range(1) depth_range(2)]);
    view(3); % Set the view to 3D
    grid on;
    hold off;
    savefig(fig1, sprintf('%s.fig', titleText));

    % Plot the 3D scatter plot of sensor positions after optimization
    fig2 = figure;
    hold on;
    scatter3(Positions_After(:,1), Positions_After(:,2), Positions_After(:,3), 'blue', 'filled');
    scatter3(centerPosition_SOA(:,1), centerPosition_SOA(:,2), centerPosition_SOA(:,3), 'red', 'filled');
    titleText = ['Day ', num2str(day), ' Sensor Positions After Optimization'];
    title(titleText);
    xlabel('Longitude(m)');
    ylabel('Latitude(m)');
    zlabel('Depth(m)');
    set(gca, 'ZDir', 'reverse'); % Continue with reversed Z axis for consistency
    axis([longitude_range(1) longitude_range(2) latitude_range(1) latitude_range(2) depth_range(1) depth_range(2)]);
    view(3);
    grid on;
    hold off;
    savefig(fig2, sprintf('%s.fig', titleText));
end
