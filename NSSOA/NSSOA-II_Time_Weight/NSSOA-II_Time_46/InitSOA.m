% Function:
% Initialize solution space for SOA algorithm
% 
% Input:
% depth: Current depth position information
% Depth: Depth layer values
% limit: Limit adjustment layers
% numSensors: Number of nodes
% population: Population size
% 
% Output:
% x_old: Initial solution generated by tent chaotic mapping
% depth_lb, depth_ub: Upper and lower bounds of adjustment for each node
% 
function [x_old, depth_lb, depth_ub] = InitSOA(depth, Depth, limit, numSensors, population)
    % Set upper and lower bounds of solution space
    depth_lb = zeros(1, numSensors);
    depth_ub = zeros(1, numSensors);
    for i = 1:numSensors
        % Current node position
        currentDepth = depth(i);

        % Find the current depth layer
        distances = sqrt((Depth - currentDepth).^2);
        [~, closestIdx] = min(distances);

        % Limit adjustment range
        lowerIndex = max(closestIdx - ((limit - 1) / 2), 1);  % Ensure not less than the minimum depth layer
        upperIndex = min(closestIdx + ((limit - 1) / 2), length(Depth));  % Ensure not greater than the maximum depth layer

        % Set upper and lower bounds of depth
        depth_lb(i) = Depth(lowerIndex);
        depth_ub(i) = Depth(upperIndex);
    end

    x_old = zeros(population, numSensors);
    % Tent chaotic mapping
    for j = 1:numSensors     
        % Generate chaotic values for population
        x = rand(1);
        chaoticValues = zeros(population, 1);
        for i = 1:population
            if x <= 0.5
                x = 2 * x;
            else
                x = 2 * (1 - x);
            end
            chaoticValues(i) = x;
        end
        % Map chaotic values to actual search space
        x_old(:, j) = depth_lb(j) + (depth_ub(j) - depth_lb(j)) * chaoticValues;
    end
end