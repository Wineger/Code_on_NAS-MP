% Function:
% Select a suitable solution from a set of dominating solutions, i.e., one that optimizes both objectives simultaneously
% 
% Inputs:
% fitness_F1: Fitness values for both objectives
% coverageRatio_noSOA, connectivityRatio_noSOA: Current uncovered ratio and connectivity ratio without optimization
% weight: Target weights
% 
% Outputs:
% index: Index of the selected solution
% 
function index = SelectSolution(fitness_F1, coverageRatio_noSOA, connectivityRatio_noSOA, weight)
    index = [];
    temp = -inf;  % Flag for comparison
    
    nPop = size(fitness_F1, 1);  % Population size
    fitness_F1_max = -fitness_F1;  % Convert to maximization problem

    for i = 1:nPop
        currentfitness = fitness_F1_max(i, :);  % Extract fitness values for both objectives
        if (currentfitness(1) > coverageRatio_noSOA) && (currentfitness(2) > connectivityRatio_noSOA)  % Optimize both objectives simultaneously
            % Find a suitable solution
            result = sum((currentfitness - [coverageRatio_noSOA, connectivityRatio_noSOA]) .* weight);
            if result > temp
                index = i;
                temp = result;
            end
        end
    end

    % If no solution optimizes both objectives simultaneously, use the one with the highest weighted sum
    if isempty(index)
        [~, index] = max(sum((fitness_F1_max - repmat([coverageRatio_noSOA, connectivityRatio_noSOA], nPop, 1)) .* weight, 2));
    end
end
