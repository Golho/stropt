function [X] = optimization(X_now, u, edof, k_0, V, V_0, X_max, X_min, p, alpha,N)
%OPTIMIZATION Summary of this function goes here
%   Detailed explanation goes here
nbrVars = length(X_now);
b_k = zeros(nbrVars,1);
ed = u(edof(:,2:end));
% Use density filter if N exists in function
if exist('N', 'var')
    X_filt = N*X_now;
    b_k = 1/alpha*(filterGradient(X_filt, p, k_0, u, edof, N, 1)) ...
        .* X_now .^(alpha+1); 
else
    for var=1:nbrVars
        u_j = ed(var, :)';
        b_k(var) = 1/alpha*(u_j'*p*X_now(var)^(p-1)*k_0{var}*u_j) * X_now(var)^(alpha+1);
    end
end
dphi_lambda_zero = @(lambda) dphi_dlambda(lambda, X_now,V, V_0, X_max, X_min, alpha,b_k);
lambda = fzero(dphi_lambda_zero, 1);
X = X_star(lambda, X_now, X_max, X_min, V, alpha, b_k);
end

function [d_phi] = dphi_dlambda(lambda, X_now, V, V_0, X_max, X_min, alpha ,b_k)
    X = X_star(lambda, X_now, X_max, X_min, V, alpha, b_k);
    d_phi = lambda * (V*sum(X) - V_0);
end

function [X] = X_star(lambda, X_now, X_max, X_min, V, alpha,b_k)
    nbrVars = length(X_now);
    X = zeros(nbrVars, 1);
    if lambda <= 0
        X = X_min;
        return;
    end
    for var = 1:nbrVars
        X(var) = (alpha*b_k(var) / (V*lambda))^(1/(alpha+1));
        if X(var) > X_max(var)
            X(var) = X_max(var);
        elseif X(var) < X_min(var)
            X(var) = X_min(var);
        end
    end
end
