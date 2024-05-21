%% Start Parallel Pool

tic

clear; clc
close all

delete(gcp('nocreate'));
par = parpool;

%% Import Data

loopStartTime = toc;

filename = 'OriginData.nc';
ImportData(filename);
filename = 'PredictData.nc';
ImportData(filename);

load OriginData.mat
load PredictData.mat
load SensorPositions.mat

disp('Data preprocessing and importing completed')
loopEndTime = toc;
disp(['Elapsed time: ', num2str((loopEndTime - loopStartTime)), ' seconds']);

Longitude = unique(data(:,1));  % Layer values
Latitude = unique(data(:,2));
Depth = unique(data(:,3));

longitude_range = [min(Longitude), max(Longitude)];  % Longitude range
latitude_range = [min(Latitude), max(Latitude)];     % Latitude range
depth_range = [min(Depth), max(Depth)];              % Depth range

u_real = U_real;
v_real = V_real;
u_predict = U_predict;
v_predict = V_predict;
u_current_ocean = u_real(:, 1);  % Initial ocean current velocity
v_current_ocean = v_real(:, 1);

%% Parameter Setting

Cxy = 1.07;   % Motion trajectory parameters
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

weight = [0.3, 0.7];  % Objective weights

interval = 4;                 % Optimization interval
number = size(u_predict, 2);  % Optimization times
day = number / interval;      % Optimization days

limit = inf;  % Limit adjustment layers

% Optimization result matrix
result = zeros(2, day);    % Predicted with optimization index
result1 = zeros(2, day);   % Real without optimization index
result2 = zeros(2, day);   % Real with optimization index

% Adjustment plan
plan_Sensors = zeros(numSensors, day);
plan_centerPosition = zeros(size(centerPosition, 1), day);

% Initializing sensor optimization deployment
loopStartTime = toc;

% [SensorPositions, Init_reslt(1, 1), Init_reslt(2, 1), F1] = InitSensorPositions(data, u_current_ocean, v_current_ocean, longitude_range, latitude_range, depth_range, Longitude, Latitude, Depth, numSensors, sensingRadius, communicationRadius, Cxy, delta_t, coverageRatio_threshold, connectivityRatio_threshold, centerPosition);

SensorPositions_noSOA = SensorPositions;
SensorPositions_SOA = SensorPositions;
% PlotSOA(F1, 0);  % Draw F1 front
% VisualizeSensorPositions(SensorPositions_noSOA, SensorPositions_SOA, 0, longitude_range, latitude_range, depth_range, centerPosition_noSOA, centerPosition_SOA);  % Visualize initial node distribution
disp(['Sensor initial optimization deployment completed: F1 member count = ' num2str(numel(F1)), ' Initial coverage rate = ', num2str(Init_reslt(1, 1)), ' Initial connectivity rate = ', num2str(Init_reslt(2, 1))]);

loopEndTime = toc;
disp(['Elapsed time: ', num2str((loopEndTime - loopStartTime) / 60), ' minutes']);

% Save initialization position information
save('SensorPositions.mat', 'SensorPositions', 'Init_reslt', 'F1');

%% Main Algorithm

for i = 1:day
    loopStartTime = toc;

    % Extracting real ocean current velocities
    u_ocean_real = u_real(:, ((i - 1) * interval + 1) : (i * interval));
    v_ocean_real = v_real(:, ((i - 1) * interval + 1) : (i * interval));

    % Extracting predicted ocean current velocities
    u_ocean_predict = u_predict(:, ((i - 1) * interval + 1) : (i * interval));
    v_ocean_predict = v_predict(:, ((i - 1) * interval + 1) : (i * interval));

    % No optimization iteration
    [SensorPositions_noSOA, result1(1, i), result1(2, i)] = noSOA(data, SensorPositions_noSOA, u_ocean_real, v_ocean_real, longitude_range, latitude_range, depth_range, sensingRadius, communicationRadius, Cxy, delta_t, centerPosition_noSOA);

    % Calculating SOA optimization adjustment scheme
    Temp1 = SensorPositions_SOA;
    Temp2 = centerPosition_SOA;
    [SensorPositions_SOA, result(1 ,i), result(2 ,i), F1, centerPosition_SOA] = SOA(data, SensorPositions_SOA, u_ocean_predict, v_ocean_predict, longitude_range, latitude_range, depth_range, Depth, limit, sensingRadius, communicationRadius, Cxy, delta_t, result1(1, i), result1(2, i), centerPosition_SOA, weight);
    plan_Sensors(:, i) = SensorPositions_SOA(:, 3) - Temp1(:, 3);
    plan_centerPosition(:, i) = centerPosition_SOA(: ,3) - Temp2(:, 3);

    % Running SOA optimization adjustment scheme
    SensorPositions(:, 3) = SensorPositions(:, 3) + plan_Sensors(:, i);
    centerPosition(:, 3) = centerPosition(:, 3) + plan_centerPosition(: ,i);
    [SensorPositions, result2(1, i), result2(2, i)] = noSOA(data, SensorPositions, u_ocean_real, v_ocean_real, longitude_range, latitude_range, depth_range, sensingRadius, communicationRadius, Cxy, delta_t, centerPosition);
    
    % Updating optimization starting position
    SensorPositions_SOA = SensorPositions;
    centerPosition_SOA = centerPosition;

    % PlotSOA(F1, i);  % Draw F1 front
    % VisualizeSensorPositions(SensorPositions_noSOA, SensorPositions_SOA, i, longitude_range, latitude_range, depth_range, centerPosition_noSOA, centerPosition_SOA);  % Visualize node distribution
    disp(['Day ', num2str(i), ' iteration optimization completed: F1 member count = ' num2str(numel(F1)), ' Current coverage rate difference = ', num2str(result2(1 ,i) - result1(1 ,i)), ' Current connectivity rate difference = ', num2str(result2(2 ,i) - result1(2 ,i))]);  % Print iteration optimization effect information

    loopEndTime = toc;
    disp(['Elapsed time: ', num2str((loopEndTime - loopStartTime) / 60), ' minutes']);
end

%% Print Results

% Concatenating initialization index values
result1 = [Init_reslt, result1];
result2 = [Init_reslt, result2];

% Outputting optimization results
disp('Without optimization results:');
disp(result1);
disp('With optimization results:');
disp(result2);
disp('Optimization differences:');
disp(result2 - result1)
disp('Sum of optimization differences:');
disp(sum(result2 - result1, 2));

% Saving results
save('result_37.mat', 'result1', 'result2');
load result_37.mat

% Draw comparison between before and after optimization
PlotOptimizationComparison(result1, result2, weight);

%% Close Parallel Pool

totalTime = toc;
disp(['Total elapsed time: ', num2str(totalTime / (60 * 60)), ' hours']);

delete(par);
