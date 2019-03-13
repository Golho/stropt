function [p] = superdraw2(Ex, Ey, X)
%SUPERDRAW2 Draws elements with specified density
%   Ex          Global elements x-coordinates [nbrElems x 4]
%   Ey          Global elements y-coordinates [nbrElems x 4]
%   X           Element densities [nbrElems x 1]
X = min(X, 1); % To ensure X <= 1
p = patch(Ex', Ey', 1-X, 'EdgeColor', 'none');
colormap('gray');
end

