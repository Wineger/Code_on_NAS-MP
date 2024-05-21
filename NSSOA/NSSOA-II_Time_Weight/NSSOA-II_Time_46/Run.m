%% Start parallel pool

tic

clear;clc
close all

delete(gcp('nocreate'));
par = parpool;

%% Import data

loopStartTime = toc;

filename = 'OriginData.nc';
ImportData(filename);

load Data.mat
load SensorPositions.mat

disp('Data preprocessing and import completed')
loopEndTime = toc;
disp(['Time elapsed: ', num2str((loopEndTime - loopStartTime)), ' seconds']);

Longitude = unique(data(:,1));  % Layering value
Latitude = unique(data(:,2));
Depth = unique(data(:,3));

longitude_range = [min(Longitude), max(Longitude)];  % Longitude range
latitude_range = [min(Latitude), max(Latitude)];     % Latitude range
depth_range = [min(Depth), max(Depth)];              % Depth range

u_real=U;
v_real=V;
u_current_ocean = u_real(:, 1);  % Initial ocean current velocity
v_current_ocean = v_real(:, 1);

%% Parameter settings

Cxy = 1.07;   % Trajectory parameter
delta_t = 6;  % Sampling time interval

numSensors = 50;            % Number of nodes
sensingRadius = 150;        % Sensing radius
communicationRadius = 150;  % Communication radius

% Convergence nodes
centerPosition = [longitude_range(1)+(1/3)*(longitude_range(2)-longitude_range(1)), latitude_range(1)+(1/3)*(latitude_range(2)-latitude_range(1)), depth_range(1)+(1/3)*(depth_range(2)-depth_range(1))
                  longitude_range(1)+(1/3)*(longitude_range(2)-longitude_range(1)), latitude_range(1)+(2/3)*(latitude_range(2)-latitude_range(1)), depth_range(1)+(2/3)*(depth_range(2)-depth_range(1))
                  longitude_range(1)+(2/3)*(longitude_range(2)-longitude_range(1)), latitude_range(1)+(1/3)*(latitude_range(2)-latitude_range(1)), depth_range(1)+(2/3)*(depth_range(2)-depth_range(1))
                  longitude_range(1)+(2/3)*(longitude_range(2)-longitude_range(1)), latitude_range(1)+(2/3)*(latitude_range(2)-latitude_range(1)), depth_range(1)+(1/3)*(depth_range(2)-depth_range(1))];
centerPosition_noSOA = centerPosition;
centerPosition_SOA = centerPosition;

coverageRatio_threshold = 0.90;     % Initialization index threshold
connectivityRatio_threshold = 0.90;

weight = [0.4, 0.6];  % Target weight

interval = 4;              % Optimization interval
number = size(u_real, 2);  % Number of optimizations
day = number / interval;   % Number of optimization days

limit = inf;  % Limit adjustment layers

% Optimization result matrix
result1 = zeros(2, day);   % No optimization index
result2 = zeros(2, day);   % With optimization index

% Initialize sensor optimization deployment
loopStartTime = toc;

% [SensorPositions, Init_reslt(1, 1), Init_reslt(2, 1), F1] = InitSensorPositions(data, u_current_ocean, v_current_ocean, longitude_range, latitude_range, depth_range, Longitude, Latitude, Depth, numSensors, sensingRadius, communicationRadius, Cxy, delta_t, coverageRatio_threshold, connectivityRatio_threshold, centerPosition);

SensorPositions_noSOA = SensorPositions;
SensorPositions_SOA = SensorPositions;
PlotSOA(F1, 0);  % Plot F1 frontier
VisualizeSensorPositions(SensorPositions_noSOA, SensorPositions_SOA, 0, longitude_range, latitude_range, depth_range, centerPosition_noSOA, centerPosition_SOA);  % Visualize initial node distribution
disp(['Sensor initial optimization deployment completed: F1 member count = ' num2str(numel(F1)), ' Initial coverage rate = ', num2str(Init_reslt(1, 1)), ' Initial connectivity rate = ', num2str(Init_reslt(2, 1))]);

loopEndTime = toc;
disp(['Time elapsed: ', num2str((loopEndTime - loopStartTime) / 60), ' minutes']);

% Save initial position information
save('SensorPositions.mat', 'SensorPositions', 'Init_reslt', 'F1');

%% Main algorithm

for i = 1:day
    loopStartTime = toc;

    % Extract each group of ocean currents
    u_ocean = u_real(:, ((i - 1) * interval + 1) : (i * interval));
    v_ocean = v_real(:, ((i - 1) * interval + 1) : (i * interval));
    
    % No optimization iteration
    [SensorPositions_noSOA, result1(1, i), result1(2, i)] = noSOA(data, SensorPositions_noSOA, u_ocean, v_ocean, longitude_range, latitude_range, depth_range, sensingRadius, communicationRadius, Cxy, delta_t, centerPosition_noSOA);

    % Mass optimization iteration
    [SensorPositions_SOA, result2(1 ,i), result2(2 ,i), F1, centerPosition_SOA] = SOA(data, SensorPositions_SOA, u_ocean, v_ocean, longitude_range, latitude_range, depth_range, Depth, limit, sensingRadius, communicationRadius, Cxy, delta_t, result1(1, i), result1(2, i), centerPosition_SOA, weight);
    
    % PlotSOA(F1, i);  % Plot F1 frontier
    % VisualizeSensorPositions(SensorPositions_noSOA, SensorPositions_SOA, i, longitude_range, latitude_range, depth_range, centerPosition_noSOA, centerPosition_SOA);  % Visualize node distribution
    disp(['Day ', num2str(i), ' optimization iteration completed: F1 member count = ' num2str(numel(F1)), ' Current coverage rate difference = ', num2str(result2(1 ,i) - result1(1 ,i)), ' Current connectivity rate difference = ', num2str(result2(2 ,i) - result1(2 ,i))]);  % Print iteration optimization effect information

    loopEndTime = toc;
    disp(['Time elapsed: ', num2str((loopEndTime - loopStartTime) / 60), ' minutes'])
end

%% Print results

% Concatenate initialization index values
result1 = [Init_reslt, result1];
result2 = [Init_reslt, result2];

% Output optimization results
disp('No optimization result:');
disp(result1);
disp('With optimization result:');
disp(result2);
disp('Optimization difference:');
disp(result2 - result1)
disp('Sum of optimization differences:');
disp(sum(result2 - result1, 2));

% Save results
save('result_46.mat', 'result1', 'result2');
load result_46.mat

% Plot optimization before and after comparison
PlotOptimizationComparison(result1, result2, weight);

%% Close parallel pool

totalTime = toc;
disp(['Total time elapsed: ', num2str(totalTime / (60 * 60)), ' hours']);

delete(par);
