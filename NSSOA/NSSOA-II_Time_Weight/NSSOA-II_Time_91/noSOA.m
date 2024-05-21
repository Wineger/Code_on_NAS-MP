% Function:
% Iterative computation of coverage and connectivity without optimization
% 
% Input:
% data: Permutations of longitude, latitude, and depth layer values
% SensorPositions: Current sensor position information
% u_ocean, v_ocean: Current ocean current velocities
% longitude_range, latitude_range, depth_range: Longitude, latitude, and depth ranges
% sensingRadius, communicationRadius: Sensing radius and communication radius
% Cxy: Motion trajectory parameters
% delta_t: Sampling time interval
% centerPosition: Convergence node
% 
% Output:
% SensorPositions: Sensor position information influenced only by motion trajectory
% coverageRatio, connectivityRatio: Coverage and connectivity without optimization
% 
function [SensorPositions, coverageRatio, connectivityRatio] = noSOA(data, SensorPositions, u_ocean, v_ocean, longitude_range, latitude_range, depth_range, sensingRadius, communicationRadius, Cxy, delta_t, centerPosition)
    % Compute the number of time intervals
    interval = size(u_ocean, 2);
    % Initialize sensor velocities to 0 at each adjustment
    u_current_sensor = u_ocean(:, 1) - u_ocean(:, 1);
    v_current_sensor = v_ocean(:, 1) - v_ocean(:, 1);

    coverageRatio_interval = zeros(1, interval);
    connectivityRatio_interval = zeros(1, interval);

    % Calculate interval indices
    for i = 1:interval
        [coverageRatio_interval(i), connectivityRatio_interval(i)] = CalculateMetrics(SensorPositions, sensingRadius, communicationRadius, longitude_range, latitude_range, depth_range, centerPosition);
        [SensorPositions, u_current_sensor, v_current_sensor] = SensorTrajectory(data, SensorPositions, u_ocean(:, i), v_ocean(:, i), u_current_sensor, v_current_sensor, longitude_range, latitude_range, Cxy, delta_t);
    end
    
    coverageRatio = mean(coverageRatio_interval);
    connectivityRatio = mean(connectivityRatio_interval);
end
