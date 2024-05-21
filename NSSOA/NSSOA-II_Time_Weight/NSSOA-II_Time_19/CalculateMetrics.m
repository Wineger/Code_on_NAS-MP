% Function:
% Calculate coverage ratio and connectivity ratio
% 
% Input:
% SensorPositions: Current sensor position information
% sensingRadius, communicationRadius: Sensing radius, communication radius
% longitude_range, latitude_range, depth_range: Longitude, latitude, depth ranges
% centerPosition: Convergence node 
% 
% Output:
% coverageRatio, connectivityRatio: Calculated coverage ratio and connectivity ratio
% 
function [coverageRatio, connectivityRatio] = CalculateMetrics(SensorPositions, sensingRadius, communicationRadius, longitude_range, latitude_range, depth_range, centerPosition)
    numSensors = size(SensorPositions, 1); % Calculate the number of nodes
    numSamplePoints = 10000; % Number of sample points for estimating coverage ratio

    % Generate random sample points for coverage estimation  
    longitudeSamples = rand(numSamplePoints, 1) * (longitude_range(2) - longitude_range(1)) + longitude_range(1);
    latitudeSamples = rand(numSamplePoints, 1) * (latitude_range(2) - latitude_range(1)) + latitude_range(1);
    depthSamples = rand(numSamplePoints, 1) * (depth_range(2) - depth_range(1)) + depth_range(1);
    samplePoints = [longitudeSamples, latitudeSamples, depthSamples];

    % Initialize coverage count
    coveredCount = 0;

    % Check if each sample point is covered by at least one sensor
    for i = 1:numSamplePoints
        for j = 1:numSensors
            if norm(samplePoints(i,:) - SensorPositions(j,:)) <= sensingRadius
                coveredCount = coveredCount + 1;
                break; % If the current point is covered, no need to check other sensors
            end
        end
    end

    % Calculate coverage ratio
    coverageRatio = coveredCount / numSamplePoints;


    state = inf(numSensors, 1);  % Initialize the layer-wise communication variable to infinity
    
    % Calculate communication situation for the first layer
    for i = 1:numSensors
        for j = 1:size(centerPosition, 1)
            if norm((SensorPositions(i, :) - centerPosition(j, :))) < communicationRadius
                state(i) = 1;
            end
        end
    end

    k = 1;
    while 1
        currentLayer = [];
        nextLayer = [];
        for i = 1:numSensors
            if state(i) == k
                currentLayer = [currentLayer; i];  % Extract indices of sensors in the current communication layer
            elseif state(i) == inf
                nextLayer = [nextLayer; i];  % Extract indices of sensors unable to communicate
            end
        end

        % Update communication situation for the next layer
        for p = 1:length(nextLayer)
            for q = 1:length(currentLayer)
                if norm(SensorPositions(nextLayer(p), :) - SensorPositions(currentLayer(q), :)) < communicationRadius
                    state(nextLayer(p)) = k + 1;
                end
            end
        end

        % If unable to communicate
        if ~any(ismember(state, k + 1))
            break;
        end

        k = k + 1;
    end
    % Calculate connectivity ratio
    connectivityRatio = (numSensors - sum(ismember(state, inf))) / numSensors;
end
