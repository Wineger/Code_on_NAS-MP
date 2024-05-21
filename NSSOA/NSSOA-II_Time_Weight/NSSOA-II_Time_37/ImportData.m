% 功能:
% 数据预处理
% 
% 输入:
% filename:文件名
% 
% 输出:
% 无
%
function ImportData(filename)
    %% 导入数据

    % 获取文件中所有变量的信息
    info = ncinfo(filename);
    
    % 初始化一个结构体来存储变量数据
    Data = struct();
    
    % 遍历文件中的每个变量
    for i = 1:length(info.Variables)
        var_name = info.Variables(i).Name; % 获取变量名
        var_size = [info.Variables(i).Dimensions.Length]; % 获取变量的维度大小
    
        % 读取变量数据
        var_data = ncread(filename, var_name);
    
        % 将数据存储在结构体中
        Data.(var_name) = var_data;
    end
    
    % 从结构体中提取各个变量
    depth = double(Data.depth);
    latitude = double(Data.latitude);
    longitude = double(Data.longitude);
    time = Data.time;
    uo = Data.uo;
    vo = Data.vo;
    
    %% 筛选数据

    % 定义经纬度范围
    minLongitude = longitude(1);  % 最小经度
    maxLongitude = longitude(2);  % 最大经度
    minLatitude = latitude(1);    % 最小纬度
    maxLatitude = latitude(2);  % 最大纬度
    
    % 找到符合条件的经纬度索引
    validLonIndices = find(longitude >= minLongitude & longitude <= maxLongitude);
    validLatIndices = find(latitude >= minLatitude & latitude <= maxLatitude);
    
    % 直接更新结构体中的经纬度数据
    Data.longitude = longitude(validLonIndices);
    Data.latitude = latitude(validLatIndices);
    
    % 更新相关的uo和vo数据
    Data.uo = Data.uo(validLonIndices, validLatIndices, :, :);
    Data.vo = Data.vo(validLonIndices, validLatIndices, :, :);
    
    %% 插值

    % 假设这是目标插值网格
    targetLon = linspace(min(Data.longitude(:)), max(Data.longitude(:)), 50);
    targetLat = linspace(min(Data.latitude(:)), max(Data.latitude(:)), 50);
    [gridLon, gridLat] = meshgrid(targetLon, targetLat);
    
    % 获取数据的维度
    [numLon, numLat, numDepth, numTime] = size(Data.uo);
    
    % 初始化插值后的结果数组
    uo_interp = zeros(length(targetLon), length(targetLat), numDepth, numTime);
    vo_interp = zeros(length(targetLon), length(targetLat), numDepth, numTime);
    
    for t = 1:numTime
        for d = 1:numDepth
            % 插值每个深度和时间的数据
            uo_interp(:, :, d, t) = interp2(Data.longitude, Data.latitude, squeeze(Data.uo(:, :, d, t)), gridLon, gridLat, 'linear');
            vo_interp(:, :, d, t) = interp2(Data.longitude, Data.latitude, squeeze(Data.vo(:, :, d, t)), gridLon, gridLat, 'linear');
        end
    end
    
    %% 更新数据

    % 更新结构体中的经纬度数据为插值网格的经纬度
    Data.longitude = gridLon(1, :)';  % 获取网格的每一列的第一个纬度值
    Data.latitude = gridLat(:, 1);   % 获取网格的每一行的第一个经度值
    
    % 更新结构体中的uo和vo数据为插值后的结果
    Data.uo = uo_interp;
    Data.vo = vo_interp;
    
    %% 拼接longitude、latitude、depth

    % 提取变量
    longitude = double(Data.longitude);
    latitude = double(Data.latitude);
    depth = double(Data.depth);
    
    % 使用 ndgrid 生成所有可能的组合
    [Lon, Lat, Dep] = ndgrid(longitude, latitude, depth);
    
    % 将生成的网格转换为列向量，这样每个点的所有组合都在一行
    Lon = Lon(:);
    Lat = Lat(:);
    Dep = Dep(:);
    
    % 组合这三个列向量成一个矩阵
    data = [Lon, Lat, Dep];
    
    %% 提取uo和vo数据

    % 确定时间的长度
    numTimes = length(time);
    
    % 初始化U和V数组
    U = zeros(size(data, 1), numTimes);
    V = zeros(size(data, 1), numTimes);
    
    % 循环遍历data中的每个位置
    for idx = 1:size(data, 1)
        % 提取当前位置的经度、纬度和深度
        current_lon = data(idx, 1);
        current_lat = data(idx, 2);
        current_depth = data(idx, 3);
        
        % 找到longitude, latitude, depth中的匹配索引
        lon_index = find(abs(longitude - current_lon) < 1e-5, 1); % 使用精度阈值防止浮点误差
        lat_index = find(abs(latitude - current_lat) < 1e-5, 1);
        depth_index = find(abs(depth - current_depth) < 1e-5, 1);
        
        % 验证索引并提取数据
        if ~isempty(lon_index) && ~isempty(lat_index) && ~isempty(depth_index)
            % 索引到uo和vo数组提取对应的所有时间点的数据
            U(idx, :) = squeeze(Data.uo(lon_index, lat_index, depth_index, :));
            V(idx, :) = squeeze(Data.vo(lon_index, lat_index, depth_index, :));
        end
    end
    
    %% 转化单位

    Longitude = unique(data(:,1));  % 分层值
    Latitude = unique(data(:,2));
    Depth = unique(data(:,3));
    
    longitude_range = [min(Longitude), max(Longitude)];  % 经度范围
    latitude_range = [min(Latitude), max(Latitude)];     % 纬度范围
    depth_range = [min(Depth), max(Depth)];              % 深度范围
    
    data(:,1) = data(:,1) * ((9*1000)/(longitude_range(2)-longitude_range(1))) * (600/(9*1000));  % 放缩区域
    data(:,2) = data(:,2) * ((9*1000)/(latitude_range(2)-latitude_range(1))) * (600/(9*1000));
    data(:,3) = data(:,3);
    
    U = U * (600/(9*1000));  % 放缩速度
    V = V * (600/(9*1000));
    
    %% 导出数据
    
    save('Data.mat', 'data', 'U', 'V');

end
