function [N] = findNeighbors(Ex, Ey, R)
%FINDNEIGHBORS Finds element neighbors and weight them
%   Ex      Global element x-coordinates
%   Ey      Global element y-coordinates
%   R       Radius to define if a element is a neighbor to another
%
%   N       [nbrElems x nbrElems] Matrix with the weighting w_ij between
%           elements i och j. The sum of every row = 1.

W = diff(Ex(1, [1 2]));
H = diff(Ey(1, [2 3]));
maxNeighbors = 1 + 4*R^2/(W*H);
nbrElems = length(Ex);
N = sparse(1:nbrElems, 1:nbrElems, R*ones(1, nbrElems), ...
    nbrElems, nbrElems, 2*nbrElems*ceil(maxNeighbors));
for current = 1:nbrElems
    xCurrent = sum(Ex(current, :))/4;
    yCurrent = sum(Ey(current, :))/4;
    for other = current+1:nbrElems
        xOther = sum(Ex(other, :))/4;
        yOther = sum(Ey(other, :))/4;
        dist = sqrt((xCurrent-xOther)^2 + (yCurrent - yOther)^2);
        weight = R - dist;
        if weight > 0
            N(current, other) = weight;
            N(other, current) = weight;
        end
    end
end
[i, j, s] = find(N);
for current = 1:nbrElems
    index = i == current;
    s(index) = s(index)/sum(s(index));
end
N = sparse(i, j, s, nbrElems, nbrElems);
end

