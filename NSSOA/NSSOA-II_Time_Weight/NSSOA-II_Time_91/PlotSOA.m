% Function:
% Plot the F1 Pareto front
%
% Inputs:
% pop: Population
% day: Current optimization day
%
% Outputs:
% None
%
function PlotSOA(pop, day)

    fig = figure;

    % Extract cost data
    Costs = -[pop.Cost];
    
    % Plot a scatter plot
    plot(Costs(1, :), Costs(2, :), 'r*', 'MarkerSize', 8);
    
    % Set the figure title and axis labels
    titleText = ['Day ', num2str(day), ' F1 Pareto front'];
    title(titleText);
    xlabel('Coverage Ratio');
    ylabel('Connectivity Ratio');
    grid on;
    
    savefig(fig, sprintf('%s.fig', titleText));
end
