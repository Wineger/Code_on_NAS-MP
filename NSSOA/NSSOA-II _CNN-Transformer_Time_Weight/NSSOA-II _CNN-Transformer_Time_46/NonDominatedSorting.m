% Function:
% Calculate non-dominated sorting of a population
% 
% Input:
% pop: Population data
% 
% Output:
% pop: Population data updated with non-dominated sorting
% F: Front data
%
function [pop, F] = NonDominatedSorting(pop)

    nPop = numel(pop);  % Number of individuals in the population

    % Initialize DominationSet and DominatedCount for each individual
    for i = 1:nPop
        pop(i).DominationSet = [];  % Indices of individuals dominated by this individual
        pop(i).DominatedCount = 0;  % Count of individuals dominating this individual
    end
    
    F{1} = [];  % Initialize Front 1 (non-dominated individuals)
    
    % Perform dominance comparison between each pair of individuals
    for i = 1:nPop
        for j = i+1:nPop
            p = pop(i);  % Individual i
            q = pop(j);  % Individual j
            
            % If p dominates q, update domination set and dominated count
            if Dominates(p, q)
                p.DominationSet = [p.DominationSet j];  % Add j to domination set of p
                q.DominatedCount = q.DominatedCount + 1;  % Increment dominated count of q
            end
            
            % If q dominates p, update domination set and dominated count
            if Dominates(q.Cost, p.Cost)
                q.DominationSet = [q.DominationSet i];  % Add i to domination set of q
                p.DominatedCount = p.DominatedCount + 1;  % Increment dominated count of p
            end
            
            pop(i) = p;  % Update individual i in the population
            pop(j) = q;  % Update individual j in the population
        end
        
        % If individual i is non-dominated, add it to Front 1 and assign rank 1
        if pop(i).DominatedCount == 0
            F{1} = [F{1} i];  % Add index of individual i to Front 1
            pop(i).Rank = 1;  % Assign rank 1 to individual i
        end
    end
    
    k = 1;  % Initialize Front index
    
    % Perform multi-level non-dominated sorting
    while true
        
        Q = [];  % Initialize empty set Q for next Front
        
        % Explore individuals in Front k
        for i = F{k}
            p = pop(i);  % Individual i
            
            % Explore individuals dominated by individual i
            for j = p.DominationSet
                q = pop(j);  % Individual j
                
                q.DominatedCount = q.DominatedCount - 1;  % Decrement dominated count of individual j
                
                % If individual j becomes non-dominated, add it to Q and assign rank k+1
                if q.DominatedCount == 0
                    Q = [Q j];  % Add index of individual j to set Q
                    q.Rank = k + 1;  % Assign rank k+1 to individual j
                end
                
                pop(j) = q;  % Update individual j in the population
            end
        end
        
        % If Q is empty, break the loop
        if isempty(Q)
            break;
        end
        
        F{k + 1} = Q;  % Assign set Q to Front k+1
        k = k + 1;  % Increment Front index
        
    end
    
end
