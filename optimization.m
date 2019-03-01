function [X] = optimization(X_now, u, edof, k_0, V, V_0, X_max, X_min, p, alpha)
%OPTIMIZATION Summary of this function goes here
%   Detailed explanation goes here
    dphi_lambda_zero = @(lambda) dphi_dlambda(lambda, X_now, u, edof, ...
        k_0, V, V_0, X_max, X_min, p, alpha);
    lambda = fzero(dphi_lambda_zero, 1);
    X = X_star(lambda, X_now, u, k_0, edof, X_max, X_min, V, p, alpha);
end

function [d_phi] = dphi_dlambda(lambda, X_now, u, edof, k_0, V, V_0, X_max, X_min, p, alpha)
    X = X_star(lambda, X_now, u, k_0, edof, X_max, X_min, V, p, alpha);
    d_phi = lambda * (V*sum(X) - V_0);
end

function [X] = X_star(lambda, X_now, u, k_0, edof, X_max, X_min, V, p, alpha)
    nbrVars = length(X_now);
    ed = u(edof(:, 2:end));
    X = zeros(nbrVars, 1);
    if lambda <= 0
        X = X_min;
        return;
    end
    for var = 1:nbrVars
        u_j = ed(var, :)';
        b_k = 1/alpha*u_j'*p*X_now(var)^(p-1)*k_0{var}*u_j * X_now(var)^(alpha+1);
        X(var) = (alpha*b_k / (V*lambda))^(1/(alpha+1));
        if X(var) > X_max(var)
            X(var) = X_max(var);
        elseif X(var) < X_min(var)
            X(var) = X_min(var);
        end
    end
end
