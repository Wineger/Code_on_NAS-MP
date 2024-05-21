% Function:
% Predict sensor trajectories
% 
% Inputs:
% data: Combination of longitude, latitude, and depth layers
% SensorPositions: Current sensor position information
% u_current_ocean, v_current_ocean: Current ocean current velocity
% u_current_sensor, v_current_sensor: Current sensor velocity
% longitude_range, latitude_range: Longitude and latitude ranges
% Cxy: Trajectory parameter
% delta_t: Sampling time interval
% 
% Outputs:
% PredictSensorPositions: Predicted sensor positions after trajectory prediction
% u_predict_sensor, v_predict_sensor: Predicted sensor velocities after trajectory prediction
% 
function [PredictSensorPositions, u_predict_sensor, v_predict_sensor] = SensorTrajectory(data, SensorPositions, u_current_ocean, v_current_ocean, u_current_sensor, v_current_sensor, longitude_range, latitude_range, Cxy, delta_t)
    % Calculate the number of nodes
    numSensors = size(SensorPositions, 1);
    
    % Initialize the predicted trajectory array to store the predicted position of each node
    PredictSensorPositions = zeros(numSensors, 3);

    % Iterate over each node for prediction
    for i = 1:numSensors
        % Initial position of the current node
        currentLongitude = SensorPositions(i, 1);
        currentLatitude = SensorPositions(i, 2);
        currentDepth = SensorPositions(i, 3);
        
        % Calculate the distance from the current position to all layers in data
        distances = sqrt((data(:,1) - currentLongitude).^2 + (data(:,2) - currentLatitude).^2 + (data(:,3) - currentDepth).^2);
        
        % Find the index of the closest node
        [~, closestIdx] = min(distances);

        % Update sensor position
        predictLongitude = currentLongitude + u_current_sensor(closestIdx) * delta_t + 1 / 2 * ((u_current_ocean(closestIdx) - u_current_sensor(closestIdx)) * Cxy) * (delta_t ^ 2);
        predictLatitude  = currentLatitude  + v_current_sensor(closestIdx) * delta_t + 1 / 2 * ((v_current_ocean(closestIdx) - v_current_sensor(closestIdx)) * Cxy) * (delta_t ^ 2);
        predictDepth = currentDepth;

        % Boundary constraint
        if predictLongitude > longitude_range(2)
            predictLongitude = longitude_range(2);
        elseif predictLongitude < longitude_range(1)
            predictLongitude = longitude_range(1);
        end
        if predictLatitude > latitude_range(2)
            predictLatitude = latitude_range(2);
        elseif predictLatitude < latitude_range(1)
            predictLatitude = latitude_range(1);
        end

        % Merge the prediction matrix
        PredictSensorPositions(i, :) = [predictLongitude, predictLatitude, predictDepth];

        % Update sensor velocity
        u_predict_sensor = u_current_sensor + ((u_current_ocean - u_current_sensor) * Cxy) * delta_t;
        v_predict_sensor = v_current_sensor + ((v_current_ocean - v_current_sensor) * Cxy) * delta_t;
    end
end
