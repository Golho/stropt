function [X_next] = optimization(X_now, u, Edof, k_0, V_e, V_0, X_max, X_min, p, alpha,N)
%OPTIMIZATION Solve a CONLIN approximation problem with Lagrangian duality
%   X_now       Densities at current point [nbrElems x 1]
%   u           Displacements in current point [nbrDofs x 1]
%   Edof        Coupling between elements and degrees of freedom 
%               [nbrElems x (1 + nbrElemsDofs)]
%   k_0         Cell array with element k_0 matrices {nbrElems x [8 x 8]}
%   V_e         Element volume
%   V_max       Maximal volume of the whole body
%   X_max       Maximum values for the densities [nbrElems x 1]
%   X_min       Minimum values for the densities [nbrElems x 1]
%   p           Penalization parameter
%   alpha       Exponent for the OC intermediate variable
%   N           Weighting matrix for the density filter 
%               [nbrElems x nbrElems]
%
%   X_next      The solution to the CONLIN approx. problem [nbrElems x 1]
X_filt = N*X_now;
b_k = 1/alpha*(filterGradient(X_filt, u, u, p, k_0, Edof, N)) ...
    .* X_now .^(alpha+1);

dphi_lambda_zero = @(lambda) dphi_dlambda(lambda, X_now,V_e, V_0, X_max, X_min, alpha, b_k);
lambda = fzero(dphi_lambda_zero, 1);

X_next = X_star(lambda, X_now, X_max, X_min, V_e, alpha, b_k);
end

function [d_phi] = dphi_dlambda(lambda, X_now, V, V_0, X_max, X_min, alpha ,b_k)
    X = X_star(lambda, X_now, X_max, X_min, V, alpha, b_k);
    d_phi = lambda * (V*sum(X) - V_0);
end

function [X] = X_star(lambda, X_now, X_max, X_min, V, alpha,b_k)
    nbrElems = length(X_now);
    X = zeros(nbrElems, 1);
    if lambda <= 0
        X = X_min;
        return;
    end
    for var = 1:nbrElems
        X(var) = (alpha*b_k(var) / (V*lambda))^(1/(alpha+1));
        if X(var) > X_max(var)
            X(var) = X_max(var);
        elseif X(var) < X_min(var)
            X(var) = X_min(var);
        end
    end
end
