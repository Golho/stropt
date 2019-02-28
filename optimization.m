function [A] = optimization(A_now, u, edof, k_0, l, V_0, A_max)
%OPTIMIZATION Summary of this function goes here
%   Detailed explanation goes here
    dphi_lambda_zero = @(lambda) dphi_dlambda(lambda, A_now, u, edof, ...
        k_0, l, V_0, A_max);
    lambda = fzero(dphi_lambda_zero, 1);
    A = A_star(lambda, A_now, u, k_0, edof, A_max, l);
end

function [d_phi] = dphi_dlambda(lambda, A_now, u, edof, k_0, l, V_0, A_max)
    A = A_star(lambda, A_now, u, k_0, edof, A_max, l);
    d_phi = lambda * (l'*A - V_0);
end

function [A] = A_star(lambda, A_now, u, k_0, edof, A_max, l)
    nbrVars = length(A_now);
    ed = u(edof(:, 2:end));
    A = 1e-4*ones(nbrVars, 1);
    if lambda > 0
       for var = 1:nbrVars
            u_j = ed(var, :)';
            A(var) = sqrt(u_j'*k_0{var}*u_j*A_now(var)^2/(l(var)*lambda));
            if A(var) > A_max(var)
                A(var) = A_max(var);
            elseif A(var) < 1e-4
                A(var) = 1e-4;
            end
        end 
    end
end
