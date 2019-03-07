function [gVal] = numGrad(g, X, h, Edof, ndof, p, k_0, F, bc, I, k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here    
nel = length(X);
gVal = zeros(nel, 1);

for e = 1:nel
    % OBS; INGEN FILTRERING SKER
    X_filt = X;
    X_filt(e) = X(e) + h;
    K = sysK(k_0, X_filt.^p, Edof, ndof);
    K = sparse(K);
    K(I ~= 0) = K(I ~= 0) + k;
    u = solveq(K, F, bc);
    gVal(e, :) = g(u);
end
end

