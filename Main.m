clear all 
clc
%%
E = 3e9;
v = 0.4;
t = 0.01;
ptype = 1; % Plane stress
w = 0.3;
h = 0.1;
nx = 3*30;
ny = 3*10;
V = w*h/(nx*ny);
V_0 = 0.4*w*h;
[Edof,Enod,coord,Dof,F,I,bc] = genMesh(w, h, nx, ny);
nel = length(Enod(:,1));
ndof = numel(Dof);

X = 0.4*ones(nel, 1);
X_max = 1*ones(nel, 1);
X_min = 1e-3*ones(nel, 1);

[Ex, Ey] = coordxtr(Edof,coord,Dof,4);
figure
eldraw2(Ex,Ey,[2,4,1]);

tol = 1e-4;
D = [E/(1-v^2),     v*E/(1-v^2),    0
    v*E/(1-v^2),    E/(1-v^2),      0
    0,              0,              E*(1-v)/(2*(1-v^2))];
%% objective function
g_0 = @(u) F'*u;
g_0plot = [];
%% k_0
k_0 = cell(nel, 1);
for e = 1:nel
    ex = Ex(e, :);
    ey = Ey(e, :);
    k_0{e} = planre(ex([1, 3]), ey([1, 3]), [ptype t], D);
end
%% Optimization loop
c = 0; % counter
X_norm = 1;
alpha = 1.5;
F = sparse(F);
while X_norm > tol
    p = 1;
    c = c + 1;
    K = sysK(k_0, X.^p, Edof, ndof);
    tic
    K = sparse(K);
    u = solveq(K, F, bc);
    toc
    g_0plot = [g_0plot g_0(u)];
    X_next = optimization(X, u, Edof, k_0, V, V_0, X_max, X_min, p, alpha);

    X_norm = norm(X_next - X);
    X = X_next;
    disp(c);
end
figure
superdraw2(Ex, Ey, X);
X_norm = 1;
while X_norm > tol && c < 600
    p = 3;
    c = c + 1;
    K = sysK(k_0, X.^p, Edof, ndof);
    u = solveq(K, F, bc);
    g_0plot = [g_0plot g_0(u)];
    X_next = optimization(X, u, Edof, k_0, V, V_0, X_max, X_min, p, alpha);

    X_norm = norm(X_next - X);
    X = X_next;
    disp(c);
end
%% Post
figure;
superdraw2(Ex, Ey, X);
figure;
plot(g_0plot);