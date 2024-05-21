%% Import Data

clear;clc
close all

filename = 'OriginData.nc';
info = ncinfo(filename);
Data = struct();

% Iterate through each variable in the file
for i = 1:length(info.Variables)
    var_name = info.Variables(i).Name;
    var_size = [info.Variables(i).Dimensions.Length];
    var_data = ncread(filename, var_name);
    Data.(var_name) = var_data;
end

% Extract variables from the structure
depth = double(Data.depth);
latitude = double(Data.latitude);
longitude = double(Data.longitude);
time = Data.time;
uo = Data.uo;
vo = Data.vo;

%% Combine longitude, latitude, and depth

% Extract variables
longitude = double(Data.longitude);
latitude = double(Data.latitude);
depth = double(Data.depth);

% Use ndgrid to generate all possible combinations
[Lon, Lat, Dep] = ndgrid(longitude, latitude, depth);

Lon = Lon(:);
Lat = Lat(:);
Dep = Dep(:);
data = [Lon, Lat, Dep];

%% Extract uo and vo data

numTimes = length(time);

U = zeros(size(data, 1), numTimes);
V = zeros(size(data, 1), numTimes);

% Loop through each position in data
for idx = 1:size(data, 1)
    % Extract the current longitude, latitude, and depth
    current_lon = data(idx, 1);
    current_lat = data(idx, 2);
    current_depth = data(idx, 3);
    
    % Find matching indices in longitude, latitude, and depth
    lon_index = find(abs(longitude - current_lon) < 1e-5, 1);
    lat_index = find(abs(latitude - current_lat) < 1e-5, 1);
    depth_index = find(abs(depth - current_depth) < 1e-5, 1);
    
    % Index into uo and vo arrays to extract data for all time points
    U(idx, :) = squeeze(uo(lon_index, lat_index, depth_index, :));
    V(idx, :) = squeeze(vo(lon_index, lat_index, depth_index, :));
end

%% Plot velocity variation with depth

flag1 = inf;
flag2 = inf;
temp1 = [];
temp2 = [];
temp3 = [];

for i = 1:length(longitude)
    for j = 1:length(latitude)
        for k = 1:length(time)
            % Extract the longitude and latitude combination
            lon = longitude(i);
            lat = latitude(j);
            tim = k;
            
            % Find all indices matching this longitude and latitude combination
            indices = find(data(:, 1) == lon & data(:, 2) == lat);
            
            % Extract depths corresponding to these indices and sort them by depth
            [sorted_depths, sorted_idx] = sort(data(indices, 3));
            
            % Use sorted indices to extract data from U and V for the first time point
            sorted_U = U(indices(sorted_idx), tim);
            sorted_V = V(indices(sorted_idx), tim);

            if (abs(sorted_U(end)) < flag1) && (abs(sorted_V(end)) < flag2)
                flag1 = sorted_U(end);
                flag2 = sorted_V(end);

                temp1 = i;
                temp2 = j;
                temp3 = k;
            end
        end
    end
end

% Extract the longitude and latitude combination
lon = longitude(temp1);
lat = latitude(temp2);
tim = temp3;

save("location.mat",'lon','lat','tim');
load location.mat

% Find all indices matching this longitude and latitude combination
indices = find(data(:, 1) == lon & data(:, 2) == lat);

% Extract depths corresponding to these indices and sort them by depth
[sorted_depths, sorted_idx] = sort(data(indices, 3));

% Use sorted indices to extract data from U and V for the first time point
sorted_U = U(indices(sorted_idx), tim);
sorted_V = V(indices(sorted_idx), tim);

% Plot velocity variation with depth
PlotVelocityDepth(sorted_depths, -sorted_U, sorted_V);
