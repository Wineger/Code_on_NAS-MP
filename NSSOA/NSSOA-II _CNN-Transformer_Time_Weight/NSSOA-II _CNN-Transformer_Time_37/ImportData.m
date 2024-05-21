% Function:
% Data preprocessing
% 
% Input:
% filename: Name of the file containing data
% 
% Output:
% None
%
function ImportData(filename)
    %% Import data

    % Get information about all variables in the file
    info = ncinfo(filename);
    
    % Initialize a structure to store variable data
    Data = struct();
    
    % Iterate over each variable in the file
    for i = 1:length(info.Variables)
        var_name = info.Variables(i).Name; % Get variable name
        var_size = [info.Variables(i).Dimensions.Length]; % Get variable dimensions
    
        % Read variable data
        var_data = ncread(filename, var_name);
    
        % Store the data in the structure
        Data.(var_name) = var_data;
    end
    
    % Extract various variables from the structure
    depth = double(Data.depth);
    latitude = double(Data.latitude);
    longitude = double(Data.longitude);
    time = Data.time;
    uo = Data.uo;
    vo = Data.vo;
    
    %% Data filtering

    % Define longitude and latitude ranges
    minLongitude = longitude(1);  % Minimum longitude
    maxLongitude = longitude(2);  % Maximum longitude
    minLatitude = latitude(1);    % Minimum latitude
    maxLatitude = latitude(2);  % Maximum latitude
    
    % Find indices of longitude and latitude that meet the conditions
    validLonIndices = find(longitude >= minLongitude & longitude <= maxLongitude);
    validLatIndices = find(latitude >= minLatitude & latitude <= maxLatitude);
    
    % Update longitude and latitude data in the structure
    Data.longitude = longitude(validLonIndices);
    Data.latitude = latitude(validLatIndices);
    
    % Update related uo and vo data
    Data.uo = Data.uo(validLonIndices, validLatIndices, :, :);
    Data.vo = Data.vo(validLonIndices, validLatIndices, :, :);
    
    %% Interpolation

    % Suppose this is the target interpolation grid
    targetLon = linspace(min(Data.longitude(:)), max(Data.longitude(:)), 50);
    targetLat = linspace(min(Data.latitude(:)), max(Data.latitude(:)), 50);
    [gridLon, gridLat] = meshgrid(targetLon, targetLat);
    
    % Get the dimensions of the data
    [numLon, numLat, numDepth, numTime] = size(Data.uo);
    
    % Initialize interpolated result arrays
    uo_interp = zeros(length(targetLon), length(targetLat), numDepth, numTime);
    vo_interp = zeros(length(targetLon), length(targetLat), numDepth, numTime);
    
    for t = 1:numTime
        for d = 1:numDepth
            % Interpolate data for each depth and time
            uo_interp(:, :, d, t) = interp2(Data.longitude, Data.latitude, squeeze(Data.uo(:, :, d, t)), gridLon, gridLat, 'linear');
            vo_interp(:, :, d, t) = interp2(Data.longitude, Data.latitude, squeeze(Data.vo(:, :, d, t)), gridLon, gridLat, 'linear');
        end
    end
    
    %% Update data

    % Update longitude and latitude data in the structure to match the interpolation grid
    Data.longitude = gridLon(1, :)';  % Get the first latitude value of each column of the grid
    Data.latitude = gridLat(:, 1);   % Get the first longitude value of each row of the grid
    
    % Update uo and vo data in the structure to the interpolated results
    Data.uo = uo_interp;
    Data.vo = vo_interp;
    
    %% Concatenate longitude, latitude, depth

    % Extract variables
    longitude = double(Data.longitude);
    latitude = double(Data.latitude);
    depth = double(Data.depth);
    
    % Use ndgrid to generate all possible combinations
    [Lon, Lat, Dep] = ndgrid(longitude, latitude, depth);
    
    % Concatenate these three column vectors into one matrix
    data = [Lon(:), Lat(:), Dep(:)];
    
    %% Extract uo and vo data

    % Determine the length of time
    numTimes = length(time);
    
    % Initialize U and V arrays
    U = zeros(size(data, 1), numTimes);
    V = zeros(size(data, 1), numTimes);
    
    % Loop through each position in data
    for idx = 1:size(data, 1)
        % Extract the longitude, latitude, and depth of the current position
        current_lon = data(idx, 1);
        current_lat = data(idx, 2);
        current_depth = data(idx, 3);
        
        % Find the matching indices in longitude, latitude, and depth
        lon_index = find(abs(longitude - current_lon) < 1e-5, 1); % Use precision threshold to prevent floating-point errors
        lat_index = find(abs(latitude - current_lat) < 1e-5, 1);
        depth_index = find(abs(depth - current_depth) < 1e-5, 1);
        
        % Validate indices and extract data
        if ~isempty(lon_index) && ~isempty(lat_index) && ~isempty(depth_index)
            % Index into uo and vo arrays to extract data for all time points
            U(idx, :) = squeeze(Data.uo(lon_index, lat_index, depth_index, :));
            V(idx, :) = squeeze(Data.vo(lon_index, lat_index, depth_index, :));
        end
    end
    
    %% Unit conversion

    Longitude = unique(data(:,1));  % Layer value
    Latitude = unique(data(:,2));
    Depth = unique(data(:,3));
    
    longitude_range = [min(Longitude), max(Longitude)];  % Longitude range
    latitude_range = [min(Latitude), max(Latitude)];     % Latitude range
    depth_range = [min(Depth), max(Depth)];              % Depth range
    
    data(:,1) = data(:,1) * ((9*1000)/(longitude_range(2)-longitude_range(1))) * (600/(9*1000));  % Rescaling area
    data(:,2) = data(:,2) * ((9*1000)/(latitude_range(2)-latitude_range(1))) * (600/(9*1000));
    data(:,3) = data(:,3);
    
    U = U * (600/(9*1000));  % Rescaling speed
    V = V * (600/(9*1000));
    
    %% Export data
    if strcmp(filename, 'OriginData.nc')
        save('OriginData.mat', 'data', 'U', 'V');
    elseif strcmp(filename, 'PredictData.nc')
        save('PredictData.mat', 'data', 'U', 'V');
    end

end
