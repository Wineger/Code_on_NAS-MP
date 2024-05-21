% Function:
% Initialize node distribution using NSSOA-II
% 
% Input:
% data: Combination of longitude, latitude, and depth layer values
% u_current_ocean, v_current_ocean: Current ocean flow velocity
% longitude_range, latitude_range, depth_range: Range of longitude, latitude, and depth
% Longitude, latitude, Depth: Longitude, latitude, and depth layer values
% numSensors: Number of nodes
% sensingRadius, communicationRadius: Sensing radius, communication radius
% Cxy: Motion trajectory parameters
% delta_t: Sampling time interval
% coverageRatio_threshold, connectivityRatio_threshold: Thresholds for initial coverage ratio and connectivity ratio
% centerPosition: Convergence node
% 
% Output:
% SensorPositions: Initial sensor position information optimized by NSSOA
% coverageRatio, connectivityRatio: Initial coverage ratio and connectivity ratio optimized by NSSOA
% F1: Dominant solution
% 
function [SensorPositions, coverageRatio, connectivityRatio, F1] = InitSensorPositions(data, u_current_ocean, v_current_ocean, longitude_range, latitude_range, depth_range, Longitude, Latitude, Depth, numSensors, sensingRadius, communicationRadius, Cxy, delta_t, coverageRatio_threshold, connectivityRatio_threshold, centerPosition)
    % Initialize random distribution of sensors
    longitude = depth_range(1) + (depth_range(2) - depth_range(1)) * rand(numSensors, 1);
    latitude = longitude_range(1) + (longitude_range(2) - longitude_range(1)) * rand(numSensors, 1);
    depth = latitude_range(1) + (latitude_range(2) - latitude_range(1)) * rand(numSensors, 1);

    limit = inf;  % No restriction on adjustment
    
    % NSSOA parameters
    fc = 2;  % Frequency parameter
    u = 1;   % Spiral constant
    v = 1;
    nPop = 50;       % Population size
    max_iter = 500;  % Maximum number of iterations

    Ms = zeros(nPop, 3*numSensors);
    Cs = zeros(nPop, 3*numSensors);
    Ds = zeros(nPop, 3*numSensors);

    % Iteration exit parameters
    count = 0;
    count_max = 100;  % Count threshold
    
    % Generate old solutions
    [x_old, x_lb, x_ub] = InitSOA(longitude, Longitude, limit, numSensors, nPop);
    [y_old, y_lb, y_ub] = InitSOA(latitude, Latitude, limit, numSensors, nPop);
    [z_old, z_lb, z_ub] = InitSOA(depth, Depth, limit, numSensors, nPop);
    X_old = [x_old, y_old, z_old];
    X_lb = [x_lb, y_lb, z_lb];
    X_ub = [x_ub, y_ub, z_ub];

    % Calculate fitness of old solutions
    fitness_old = zeros(nPop, 2);
    for i =1:nPop
        SensorPositions_old = [X_old(i, 1 : numSensors); X_old(i, numSensors + 1 : 2 * numSensors); X_old(i, 2 * numSensors + 1 : end)]';
        fitness_old(i,:) = Obj_SOA(data,  SensorPositions_old, u_current_ocean, v_current_ocean, longitude_range, latitude_range, depth_range, sensingRadius, communicationRadius, Cxy, delta_t, centerPosition);
    end

    % Fast non-dominated sorting of old solutions
    empty_individual.Position = [];
    empty_individual.Cost = [];
    empty_individual.Rank = [];
    empty_individual.DominationSet = [];
    empty_individual.DominatedCount = [];
    empty_individual.CrowdingDistance = [];
    
    pop = repmat(empty_individual, nPop, 1);
    for i = 1:nPop
        pop(i).Position = X_old(i,:);
        pop(i).Cost = fitness_old(i,:)';
    end

    [pop, F] = NonDominatedSorting(pop);  % Non-dominated sorting    
    pop = CalcCrowdingDistance(pop, F);   % Calculate crowding distance 
    [pop, F] = SortPopulation(pop);       % Sort

    % Find F1 frontier
    F1 = pop(F{1});
    for i = 1:size(F1, 2)
        X_F1(i,:) = F1(i).Position;
        fitness_F1(i,:) = F1(i).Cost';
    end
    
    for iter = 1:max_iter

        X_F1_temp = X_F1;  % F1 frontier of the last iteration
        PbestX = X_old(1, :);  % Global best solution of the last iteration

        % Generate new solutions
        for i = 1:nPop
            % Migration behavior
            A = fc - fc * (iter / max_iter);
            Cs(i, :) =  X_old(i, :) * A;

            B = 2 * A^2 * rand(1);
            Ms(i,:) = B * (PbestX - X_old(i, :));

            Ds(i, :) = abs(Cs(i, :) + Ms(i, :));

            % Attack behavior
            theta = rand(1);
            r = u * exp(v * theta);
            x = r * cos(2 * pi * theta);
            y = r * sin(2 * pi * theta);
            z = r * theta;

            X_new(i, :) = x * y * z * Ds(i, :) + PbestX;
        end

        % Correct new solutions
        for i =1:nPop
            for j =1:3*numSensors
                if X_new(i, j) > X_ub(j)
                    X_new(i, j) = X_ub(j);
                elseif X_new(i, j) < X_lb(j)
                    X_new(i, j) = X_lb(j);
                end
            end
        end

        % Calculate fitness of new solutions
        for i = 1:nPop
            SensorPositions_new = [X_new(i, 1 : numSensors); X_new(i, numSensors + 1 : 2 * numSensors); X_new(i, 2 * numSensors + 1 : end)]';
            fitness_new(i,:) = Obj_SOA(data, SensorPositions_new, u_current_ocean, v_current_ocean, longitude_range, latitude_range, depth_range, sensingRadius, communicationRadius, Cxy, delta_t, centerPosition);
        end

        % Merge old and new solutions
        X = [X_new; X_old];
        fitness = [fitness_new; fitness_old];

        % Fast non-dominated sorting of new and old solutions
        pop = [pop; pop];
        for i = 1:2*nPop
            pop(i).Position = X(i,:);
            pop(i).Cost = fitness(i,:)';
        end
            
        [pop, F] = NonDominatedSorting(pop);  % Non-dominated sorting
        pop = CalcCrowdingDistance(pop, F);   % Calculate crowding distance 
        pop = SortPopulation(pop);            % Sort
        
        % Fast non-dominated sorting of new solutions
        pop = pop(1:nPop);  % Extract new solutions

        [pop, F] = NonDominatedSorting(pop);  % Non-dominated sorting 
        pop = CalcCrowdingDistance(pop, F);   % Calculate crowding distance      
        [pop, F] = SortPopulation(pop);       % Sort

        % Update old solutions and fitness
        for i = 1:nPop
            X_old(i,:) = pop(i).Position;
            fitness_old(i,:) = pop(i).Cost';
        end

        % Find F1 frontier
        F1 = pop(F{1});
        for i = 1:size(F1, 2)
            X_F1(i,:) = F1(i).Position;
            fitness_F1(i,:) = F1(i).Cost';
        end

        % Automatic iteration exit
        if X_F1 == X_F1_temp
            count = count + 1;
        else
            count = 0;
        end
        if count > count_max
            break;
        end
    end
   
    % Find appropriate solutions
    temp = -inf;  % Judgment flag
    weight = [0.3, 0.7];  % Target weights

    for i = 1:nPop
        currentfitness = -pop(i).Cost';  % Extract fitness of two objectives in turn
        if ((currentfitness(1) < coverageRatio_threshold) && (currentfitness(2) < connectivityRatio_threshold))  % Satisfy both thresholds

            result = sum(currentfitness .* weight);
            if result > temp
                Solution_index = i;
                temp = result;
            end

        end
    end

    coverageRatio = -pop(Solution_index).Cost(1);
    connectivityRatio = -pop(Solution_index).Cost(2);
    SensorPositions = [pop(Solution_index).Position(1 : numSensors); pop(Solution_index).Position(numSensors + 1 : 2 * numSensors); pop(Solution_index).Position(2 * numSensors + 1 : end)]';
end
