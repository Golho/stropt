function [K] = sysK(k_0, E, edof, nbrDofs)
%SYSK Calculate global stiffness matric from element k_0-submatrices
%   k_0         Cell array with element k_0 matrices {nbrElems x [8 x 8]}
%   E           E modulus for each element [nbrElems x 1]
%   edof        Coupling between elements and degrees of freedom
%   nbrDofs     Number of degrees of freedom
%
%   K           Global stiffness matrix [nbrDofs x nbrDofs]
nbrElem = length(E);
K = zeros(nbrDofs, nbrDofs);
for e = 1:nbrElem
    Ke = E(e)*k_0{e};
    index = edof(e, 2:end);
    K(index, index) = K(index, index) + Ke;
end
end

