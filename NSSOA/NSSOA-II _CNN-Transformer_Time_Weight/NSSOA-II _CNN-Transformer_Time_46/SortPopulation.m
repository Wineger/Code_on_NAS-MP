% Function:
% Sort the population in ascending order
% 
% Input:
% pop: population
% 
% Output:
% pop: sorted population
% F: front data
% 
function [pop, F] = SortPopulation(pop)

    % Sort by crowding distance
    [~, CDSO] = sort([pop.CrowdingDistance], 'descend');
    pop = pop(CDSO);
    
    % Sort by rank
    [~, RSO] = sort([pop.Rank]);
    pop = pop(RSO);
    
    % Update fronts
    Ranks = [pop.Rank];
    MaxRank = max(Ranks);
    F = cell(MaxRank, 1);
    for r = 1:MaxRank
        F{r} = find(Ranks == r);
    end

end
