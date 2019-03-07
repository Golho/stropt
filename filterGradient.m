function [gradient] = filterGradient(X_filt, u_1, u_2, p, k_0, edof, N)
% Calculate u_1' * dK/drho_j * u_2 where K = X_filt^p * k_0
%   X_filt      Filtered densities [nbrElems x 1]
%   u_1         Displacements u_1 (see above) [nbrDofs x 1]
%   u_2         Displacements u_2 (see above) [nbrDofs x 1]
%   p           Penalization parameter
%   k_0         Cell array with element k_0 matrices {nbrElems x [8 x 8]}
%   edof        Coupling between elements and degrees of freedom
%   N           Weighting matrix for the density filter 
%               [nbrElems x nbrElems]
%
%   gradient    Gradient (see above) [nbrElems x 1]
nbrElems = length(X_filt);
gradient = zeros(nbrElems, 1);
for e = 1:nbrElems
    idof = edof(e, 2:end);
    u_1j = u_1(idof);
    u_2j = u_2(idof);
    neighbors = find(N(e, :));
    for iNeighbor = neighbors
        w = N(e, iNeighbor);
        gradient(e) = gradient(e) + w * u_1j'* p*X_filt(e)^(p-1)*k_0{e} *u_2j;
    end
end
end

