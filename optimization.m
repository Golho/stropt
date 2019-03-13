function [X_next] = optimization(X_min, X_max, g_1OC, X_star_fun)
%OPTIMIZATION Solve a CONLIN approximation problem with Lagrangian duality
%   X_min       Minimum values for the densities [nbrElems x 1]
%   X_max       Maximum values for the densities [nbrElems x 1]
%   g_1OC       First approx. constraint as function handle @(X), X =
%               [nbrElems x 1]
%   X_star_fun  Function handle for the X_star formula @(lambda)
%
%   X_next      The solution to the CONLIN approx. problem [nbrElems x 1]
dphi_lambda_zero = @(lambda) dphi_dlambda(lambda, X_min, X_max, g_1OC, X_star_fun);
lambda = fzero(dphi_lambda_zero, 1);

X_next = X_star(lambda, X_min, X_max, X_star_fun);
end

function [d_phi] = dphi_dlambda(lambda, X_min, X_max, g_1OC, X_star_fun)
    X = X_star(lambda, X_min, X_max, X_star_fun);
    d_phi = lambda * g_1OC(X);
end

function [X] = X_star(lambda, X_min, X_max, X_star_fun)
    if lambda <= 0
        X = X_min;
        return;
    end
    X = X_star_fun(lambda);
    X(X > X_max) = X_max(X > X_max);
    X(X < X_min) = X_min(X < X_min);
end
