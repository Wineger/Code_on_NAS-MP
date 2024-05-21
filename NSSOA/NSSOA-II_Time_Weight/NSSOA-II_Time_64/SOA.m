% Function:
% NSSOA-II: SOA optimization algorithm based on elite strategy for fast non-dominated sorting-II
% 
% Input:
% data: combinations of longitude, latitude, and depth layers
% SensorPositions: current sensor position information
% u_ocean, v_ocean: current ocean current speed
% longitude_range, latitude_range, depth_range: longitude, latitude, and depth ranges
% Depth: depth layer values
% limit: limit the number of layers to adjust
% sensingRadius, communicationRadius: sensing radius, communication radius
% Cxy: motion trajectory parameters
% delta_t: sampling time interval
% coverageRatio_noSOA, connectivityRatio_noSOA: current coverage rate and connectivity rate without optimization
% centerPosition: convergence node
% weight: objective weight
% 
% Output:
% SensorPositions: sensor position information after SOA optimization for interval time steps
% coverageRatio, connectivityRatio: coverage rate and connectivity rate after SOA optimization
% F1: dominant solution
% 
function [SensorPositions, coverageRatio, connectivityRatio, F1, centerPosition] = SOA(data, SensorPositions, u_ocean, v_ocean, longitude_range, latitude_range, depth_range, Depth, limit, sensingRadius, communicationRadius, Cxy, delta_t, coverageRatio_noSOA, connectivityRatio_noSOA, centerPosition, weight)
    numSensors = size(SensorPositions, 1);
    numcenterPosition = size(centerPosition, 1);
    num = numSensors + numcenterPosition;  % Calculate the number of nodes

    longitude = [SensorPositions(:, 1); centerPosition(:, 1)];  % Separate position information
    latitude = [SensorPositions(:, 2); centerPosition(:, 2)];
    depth = [SensorPositions(:, 3); centerPosition(:, 3)];

    % SOA parameters
    fc = 2;  % Frequency parameter
    u = 1;   % Spiral constant
    v = 1;
    nPop = 50;       % Population size
    max_iter = 500;  % Number of iterations

    Ms = zeros(nPop, num);
    Cs = zeros(nPop, num);
    Ds = zeros(nPop, num);

    % Iteration exit parameters
    count = 0;
    count_max = 100;  % Counting threshold
    
    % Generate old solutions
    [x_old, x_lb, x_ub] = InitSOA(depth, Depth, limit, num, nPop);

    % Calculate the fitness of old solutions
    fitness_old = zeros(nPop, 2);
    for i = 1:nPop
        SensorPositions_old = [longitude(1:numSensors), latitude(1:numSensors), x_old(i, 1:numSensors)'];
        centerPosition_old = [longitude(numSensors + 1:end), latitude(numSensors + 1:end), x_old(i, numSensors + 1:end)'];
        fitness_old(i, :) = Obj_SOA(data,  SensorPositions_old, u_ocean, v_ocean, longitude_range, latitude_range, depth_range, sensingRadius, communicationRadius, Cxy, delta_t, centerPosition_old);
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
        pop(i).Position = x_old(i, :);
        pop(i).Cost = fitness_old(i, :)';
    end

    [pop, F] = NonDominatedSorting(pop);  % Non-dominated sorting    
    pop = CalcCrowdingDistance(pop, F);   % Calculate crowding distance 
    [pop, F] = SortPopulation(pop);       % Sort
    
    % Find F1 front
    F1 = pop(F{1});
    for i = 1:size(F1, 2)
        x_F1(i, :) = F1(i).Position;
        fitness_F1(i, :) = F1(i).Cost';
    end
    
    for iter = 1:max_iter

        x_F1_temp = x_F1;  % F1 front of the previous iteration
        PbestX = x_old(1, :);  % Global optimal solution of the previous iteration

        % Generate new solutions
        for i = 1:nPop
            % Migration behavior
            A = fc - fc * (iter / max_iter);
            Cs(i, :) =  x_old(i, :) * A;

            B = 2 * A^2 * rand(1);
            Ms(i, :) = B * (PbestX - x_old(i, :));

            Ds(i, :) = abs(Cs(i, :) + Ms(i, :));

            % Attack behavior
            theta = rand(1);
            r = u * exp(v * theta);
            x = r * cos(2 * pi * theta);
            y = r * sin(2 * pi * theta);
            z = r * theta;

            x_new(i, :) = x * y * z * Ds(i, :) + PbestX;
        end

        % Correct new solutions
        for i = 1:nPop
            for j = 1:num
                if x_new(i, j) > x_ub(j)
                    x_new(i, j) = x_ub(j);
                elseif x_new(i, j) < x_lb(j)
                    x_new(i, j) = x_lb(j);
                end
            end
        end

        % Calculate the fitness of new solutions
        for i = 1:nPop
            SensorPositions_new = [longitude(1:numSensors), latitude(1:numSensors), x_new(i, 1:numSensors)'];
            centerPosition_new = [longitude(numSensors + 1:end), latitude(numSensors + 1:end), x_new(i, numSensors + 1:end)'];
            fitness_new(i, :) = Obj_SOA(data, SensorPositions_new, u_ocean, v_ocean, longitude_range, latitude_range, depth_range, sensingRadius, communicationRadius, Cxy, delta_t, centerPosition_new);
        end

        % Merge old and new solutions
        X = [x_new; x_old];
        fitness = [fitness_new; fitness_old];

        % Fast non-dominated sorting of new and old solutions
        pop = [pop; pop];
        for i = 1:2 * nPop
            pop(i).Position = X(i, :);
            pop(i).Cost = fitness(i, :)';
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
            x_old(i, :) = pop(i).Position;
            fitness_old(i, :) = pop(i).Cost';
        end

        % Find F1 front
        F1 = pop(F{1});
        for i = 1:size(F1, 2)
            x_F1(i, :) = F1(i).Position;
            fitness_F1(i, :) = F1(i).Cost';
        end

        % Automatic exit iteration
        if x_F1 == x_F1_temp
            count = count + 1;
        else
            count = 0;
        end
        if count > count_max
            break;
        end
    end

    % Find F1 front
    F1 = pop(F{1});
    for i = 1:size(F1, 1)
        x_F1(i, :) = F1(i).Position;
        fitness_F1(i, :) = F1(i).Cost';
    end
   
    % Find suitable solutions
    Solution_index = SelectSolution(fitness_F1, coverageRatio_noSOA, connectivityRatio_noSOA, weight);

    % Extract two objective values
    coverageRatio = -fitness_F1(Solution_index, 1);
    connectivityRatio = -fitness_F1(Solution_index, 2);

    % Extract position information
    depth = x_F1(Solution_index, :);
    SensorPositions = [longitude(1:numSensors), latitude(1:numSensors), depth(1:numSensors)'];
    centerPosition = [longitude(numSensors + 1:end), latitude(numSensors + 1:end), depth(numSensors + 1:end)'];
    [~, SensorPositions] = Obj_SOA(data, SensorPositions, u_ocean, v_ocean, longitude_range, latitude_range, depth_range, sensingRadius, communicationRadius, Cxy, delta_t, centerPosition);
end
