clear; clc;
close all;

load result_55.mat;
day = size(result1, 2);
x = 0 : day-1;

figure(200);
hold on;

% Store handles and labels for the plots
h = [];
labels = [];

ColorType = [144, 190, 224] ./ 255;
result_Before = sum(result1 .* repmat([0.5; 0.5], 1, day));
h(end+1) = plot(x, result_Before, '--o', 'Color', ColorType, 'LineWidth', 1.6, 'MarkerSize', 7);
labels{end+1} = 'Without optimization';

for i = 1:5
    switch i
        case 1
            load result_91.mat;
            weight = [0.9; 0.1];
            ColorType = [131, 064, 038];
        case 2
            load result_82.mat;
            weight = [0.8; 0.2];
            ColorType = [130, 178, 154];
        case 3
            load result_73.mat;
            weight = [0.7; 0.3];
            ColorType = [033, 158, 188];
        case 4
            load result_64.mat;
            weight = [0.6; 0.4];
            ColorType = [255, 223, 146];
        case 5
            load result_55.mat;
            weight = [0.5; 0.5];
            ColorType = [217, 079, 051];
    end

    ColorType = ColorType ./ 255;
    
    result_After = sum(result2 .* repmat([0.5; 0.5], 1, day));
    h(end+1) = plot(x, result_After, '-*', 'Color', ColorType,'LineWidth', 1.6, 'MarkerSize', 7);
    labels{end+1} = ['Weights ', num2str(weight(1)), ' and ', num2str(weight(2)), ' with optimization'];
end

% Reverse the handles and labels for the legend
h = fliplr(h);
labels = fliplr(labels);

xlabel('Optimization cycles');
ylabel('Weighted sum of coverage and connectivity rates');
legend(h, labels, 'Location', 'best');
grid on;
hold off;

figHandle = figure(200);
fileName = 'Change_Weights.fig';
savefig(figHandle, fileName);
