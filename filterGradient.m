function [gradient] = filterGradient(X_filt, p, k_0, u, edof, N, CASE)
% Derivate [dK/drho_i * u]
nbrVars = length(X_filt);
nbrDofs = length(u);
if CASE == 1
    gradient = zeros(nbrVars, 1);
elseif CASE == 2
    gradient = spalloc(nbrDofs, nbrVars, 8*nbrVars);
end
ed = u(edof(:, 2:end));
for var = 1:nbrVars
    varDof = edof(var, 2:end);
    dfdx = zeros(nnz(N(var,:)),1);
    c = 1;
    neighbors = find(N(var, :));
    sumGrad = zeros(8, 1);
    for ie = neighbors
        dxidxe = N(var,ie);
        u_j = ed(ie, :)';
        if CASE == 1
            dfdxi = u_j'*p*X_filt(ie)^(p-1)*k_0{ie}*u_j;
            dfdx(c) = dxidxe*dfdxi;
        elseif CASE == 2
            grad_e = p*X_filt(ie)^(p-1)*k_0{ie}*u_j;
            sumGrad = sumGrad + grad_e*dxidxe;
        end
        c = c + 1;
    end
    
    if CASE == 1
        gradient(var) = sum(dfdx);
    elseif CASE == 2
        gradient(varDof, var) = gradient(varDof, var) + sumGrad;
    end
end
end

