function [gradient] = filterGradient(dEdrho, u_1, u_2, k_0, edof, N)
% Calculate u_1' * dK/drho_j * u_2 where k_e = E*k_0
%   dEdrho      Gradient E of rho (filtered) [nbrElems x 1]
%   u_1         Displacements u_1 (see above) [nbrDofs x 1]
%   u_2         Displacements u_2 (see above) [nbrDofs x 1]
%   k_0         Cell array with element k_0 matrices {nbrElems x [8 x 8]}
%   edof        Coupling between elements and degrees of freedom
%   N           Weighting matrix for the density filter 
%               [nbrElems x nbrElems]
%
%   gradient    Gradient (see above) [nbrElems x 1]
nbrElems = length(dEdrho);
dK = zeros(nbrElems, 1);
for e = 1:nbrElems
    idof = edof(e, 2:end);
    u_1j = u_1(idof);
    u_2j = u_2(idof);
    dK(e) = u_1j' * dEdrho(e) * k_0{e} * u_2j;
end
gradient = N*dK;
end

