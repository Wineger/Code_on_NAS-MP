% Function:
% Plot comparison of optimization effects
%
% Inputs:
% result1, result2: Sensor metrics before and after optimization
% weight: Target weights
%
% Outputs:
% None
%
function PlotOptimizationComparison(result1, result2, weight)
    % Extract coverage and connectivity metrics before and after optimization
    coverageRatio_After = result2(1, :);
    connectivityRatio_After = result2(2, :);
    coverageRatio_Before = result1(1, :);
    connectivityRatio_Before = result1(2, :);

    % Calculate the weighted sum
    result_After = coverageRatio_After .* 0.5 + connectivityRatio_After .* 0.5;
    result_Before = coverageRatio_Before .* 0.5 + connectivityRatio_Before .* 0.5;

    % Calculate the optimization period
    day = length(coverageRatio_Before) - 1;
    x = 0 : day;

    % Define custom RGBcolors
    red = [217, 079, 051] ./ 255;
    blue = [144, 190, 224] ./ 255;

    % Plot coverage ratio comparison
    figure(200);
    hold on;
    plot(x, coverageRatio_After, '-*', 'Color', red, 'LineWidth', 1.6, 'MarkerSize', 7, 'DisplayName', 'Coverage After Optimization');
    plot(x, coverageRatio_Before, '-o', 'Color', blue, 'LineWidth', 1.6, 'MarkerSize', 7, 'DisplayName', 'Coverage Before Optimization');    
    xlabel('Optimization Cycles');
    ylabel('Coverage Ratio');
    legend('show');
    grid on;
    hold off;

    figHandle = figure(200);
    fileName = 'Effect of Optimization Cycles on Coverage.fig';
    savefig(figHandle, fileName);
    
    % Plot connectivity ratio comparison
    figure(201);
    hold on;
    plot(x, connectivityRatio_After, '-*', 'Color', red, 'LineWidth', 1.6, 'MarkerSize', 7, 'DisplayName', 'Connectivity After Optimization');
    plot(x, connectivityRatio_Before, '-o', 'Color', blue, 'LineWidth', 1.6, 'MarkerSize', 7, 'DisplayName', 'Connectivity Before Optimization');
    xlabel('Optimization Cycles');
    ylabel('Connectivity Ratio');
    legend('show');
    grid on;
    hold off;

    figHandle = figure(201);
    fileName = 'Effect of Optimization Cycles on Connectivity.fig';
    savefig(figHandle, fileName);

    % Plot weighted sum comparison
    figure(202);
    hold on;
    plot(x, result_After, '-*', 'Color', red, 'LineWidth', 1.6, 'MarkerSize', 7, 'DisplayName', ['Weights ', num2str(weight(1)), ' and ', num2str(weight(2)), ' With Optimization']);
    plot(x, result_Before, '-o', 'Color', blue, 'LineWidth', 1.6, 'MarkerSize', 7, 'DisplayName', 'Without Optimization');    
    xlabel('Optimization Cycles');
    ylabel('Weighted Sum of Coverage and Connectivity');
    legend('show');
    grid on;
    hold off;

    figHandle = figure(202);
    fileName = 'Effect of Optimization Cycles on Weighted Sum.fig';
    savefig(figHandle, fileName);
end
