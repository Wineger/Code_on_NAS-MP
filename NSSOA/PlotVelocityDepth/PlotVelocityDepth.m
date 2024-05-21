function PlotVelocityDepth(Depth, u, v)
    % Define custom RGB colors
    red = [217, 079, 051] ./ 255;
    blue = [144, 190, 224] ./ 255;

    % Create a new figure window
    figure(200);
    
    hold on;
    plot(u, Depth, '-', 'Color' , red, 'LineWidth', 1.6);
    plot(v, Depth, '-', 'Color' , blue, 'LineWidth', 1.6);
    
    xlabel('Velocity (m/s)');
    ylabel('Depth (m)');
    legend('Westward current velocity', 'Northward current velocity');
    set(gca, 'YDir','reverse');
    axis([(min([u;v]) - 0.02) (max([u;v]) + 0.02) min(Depth) max(Depth)])
    grid on;
    hold off;

    figHandle = figure(200);
    fileName = 'Velocity Change with Depth.fig';
    savefig(figHandle, fileName);
end
