% Function:
% Calculate dominance relationship between two solutions
% 
% Input:
% x, y: Two solutions. They can be either vectors or structures with a field 'Cost' representing the solution's cost.
% 
% Output:
% b: Dominance relationship (1 if x dominates y, 0 if y dominates x)
%
function b = Dominates(x, y)
    % Check if x is a structure, if so, extract its Cost field
    if isstruct(x)
        x = x.Cost;
    end

    % Check if y is a structure, if so, extract its Cost field
    if isstruct(y)
        y = y.Cost;
    end

    % Determine dominance relationship
    % A solution x dominates y if all elements of x are less than or equal to the corresponding elements of y,
    % and at least one element of x is strictly less than the corresponding element of y
    b = all(x <= y) && any(x<y);
end
