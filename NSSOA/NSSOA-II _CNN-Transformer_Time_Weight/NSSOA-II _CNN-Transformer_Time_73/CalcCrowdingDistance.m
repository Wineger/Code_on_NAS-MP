% Function:
% Calculates the crowding distance of the population based on the given front data
% 
% Input:
% pop: Population data
% F: Front data
% 
% Output:
% pop: Population data updated with crowding distance
%
function pop = CalcCrowdingDistance(pop, F)
    % Determine the number of fronts
    nF = numel(F);
    
    % Iterate over each front
    for k = 1:nF
        % Extract costs for individuals in the current front
        Costs = [pop(F{k}).Cost];
        
        % Determine the number of objectives
        nObj = size(Costs, 1);
        
        % Determine the number of individuals in the current front
        n = numel(F{k});
        
        % Initialize crowding distance matrix
        d = zeros(n, nObj);
        
        % Calculate crowding distance for each objective
        for j = 1:nObj
            % Sort costs for the current objective
            [cj, so] = sort(Costs(j, :));
            
            % Set boundary individuals' crowding distance to infinity
            d(so(1), j) = inf;
            d(so(end), j) = inf;
            
            % Calculate crowding distance for intermediate individuals
            for i = 2:n-1
                d(so(i), j) = abs(cj(i+1)-cj(i-1))/abs(cj(1)-cj(end));
            end
        end
        
        % Sum crowding distances for each individual and assign to the population data
        for i = 1:n
            pop(F{k}(i)).CrowdingDistance = sum(d(i, :));
        end
    end
end
