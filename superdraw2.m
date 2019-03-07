function [] = superdraw2(Ex, Ey, X)
%SUPERDRAW2 Draws elements with specified density
%   Ex          Global elements x-coordinates [nbrElems x 4]
%   Ey          Global elements y-coordinates [nbrElems x 4]
%   X           Element densities [nbrElems x 1]
nel = length(Ex);
X = min(X, 1); % To ensure X <= 1
for e = 1:nel
    ex = Ex(e, :);
    ey = Ey(e, :);
    patch(ex, ey, [1-X(e), 1-X(e), 1-X(e)]);
end
end

