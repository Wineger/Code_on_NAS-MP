% Function:
% Objective function for SOA
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
% y: Fitness values for two objectives
% SensorPositions: Position information after passing through interval time steps
% 
function [y, SensorPositions] = Obj_SOA(data, SensorPositions, u_ocean, v_ocean, longitude_range, latitude_range, depth_range, sensingRadius, communicationRadius, Cxy, delta_t, centerPosition)

    % Calculate interval indices
    [SensorPositions, coverageRatio, connectivityRatio] = noSOA(data, SensorPositions, u_ocean, v_ocean, longitude_range, latitude_range, depth_range, sensingRadius, communicationRadius, Cxy, delta_t, centerPosition);
    
    % Multi-objective optimization
    y = -[coverageRatio connectivityRatio];
end
