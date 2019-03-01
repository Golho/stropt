%% Describe coordinates and couplings
coord = [0 0;
         1 0;
         2 0;
         3 0;
         4 0;
         3 1;
         2 1;
         1 1;
         0 1];
dof = [(1:2:18)' (2:2:18)'];
enod = [1 1 2;
        2 2 3;
        3 3 4;
        4 4 5;
        5 5 6;
        6 4 6;
        7 4 7;
        8 3 7;
        9 3 8;
        10 2 8;
        11 2 9;
        12 9 8;
        13 8 7;
        14 7 6];
edof = [enod(:, 1) dof(enod(:, 2), :) dof(enod(:, 3), :)];
nbrDofs = numel(dof);
nbrElem = size(edof, 1);
%% Describe loading and boundary conditions
F = zeros(nbrDofs, 1);
F(10) = -500;
 bc = [1 0;
       2 0;
       17 0;
       18 0];
E = 210e9;
A_0 = 0.1*ones(nbrElem, 1);
A_max = 0.1*ones(nbrElem, 1);
%% Calculate displacements
[Ex, Ey] = coordxtr(edof, coord, dof, 2);
A = A_0;
A_norm = 1;
tol = 1e-6;
V_0 = 1000/7800;
l = sqrt((Ex(:, 1)-Ex(:, 2)).^2+(Ey(:, 1)-Ey(:, 2)).^2);
g_0 = @(u) F'*u;
g_1 = @(A) A*l - V_0;
%% k_0
k_0 = cell(nbrElem, 1);
for e = 1:nbrElem
    ex = Ex(e, :);
    ey = Ey(e, :);
    k_0{e} = bar2e(ex, ey, [E 1]);
end
%% Optimization loop
c = 0;
while A_norm > tol
    c = c + 1;
    close all;
    K = sysK(k_0, A, edof, nbrDofs);
    u = solveq(K, F, bc);
    A_next = optimization(A, u, edof, k_0, l, V_0, A_max);
    A_norm = norm(A_next - A);
    A = A_next;
end
%% Plot
figure;
myeldraw2(Ex, Ey, [1 2 1], A, 1e3);
myeldisp2(Ex, Ey, u(edof(:, 2:end)), [1 4 1], 1e5, A, 1e3);
%% Post-calculate stress
s = zeros(nbrElem, 1);
for e = 1:nbrElem
    ex = Ex(e,:);
    ey = Ey(e,:);
    ed = u(edof(e, 2:end))';
    s(e) = bar2s(ex, ey, [E, A(e)], ed)/A(e);
end