function [] = superdraw2(Ex, Ey, X)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
nel = length(Ex);
X = min(X, 1);
for e = 1:nel
    ex = Ex(e, :);
    ey = Ey(e, :);
    patch(ex, ey, [1-X(e), 1-X(e), 1-X(e)]);
end
end

