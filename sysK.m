function [K] = sysK(k_0, X, edof, nbrDofs)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
nbrElem = length(X);
K = zeros(nbrDofs, nbrDofs);
for e = 1:nbrElem
    Ke = X(e)*k_0{e};
    index = edof(e, 2:end);
    K(index, index) = K(index, index) + Ke;
end
end

