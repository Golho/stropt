clear all; clc; close all;
%% OPTIONS
MESH_SCALE = 20;
P_PAR = 3; % p-parameter
ALPHA = 1.5; % alpha-parameter
RADIUS = 0.01; % Radius of Weight function
%% Initializing parameters and generating mesh
% Material parameters
E = 3e9;
v = 0.4;
t = 0.01;
ptype = 1; % Plane stress
D = [E/(1-v^2),     v*E/(1-v^2),    0
     v*E/(1-v^2),    E/(1-v^2),      0
     0,              0,              E*(1-v)/(2*(1-v^2))];

% Generate mesh
w = 0.3;
h = 0.1;
nx = MESH_SCALE*3;
ny = MESH_SCALE*1;
[Edof, Enod, coord, Dof, F, I, bc] = genMesh(w, h, nx, ny);
F = sparse(F);
[Ex, Ey] = coordxtr(Edof,coord,Dof,4);
%N = findNeighbors(Ex, Ey, 0.021);
nel = length(Enod(:,1));
ndof = numel(Dof);

V = w*h/(nx*ny); % element volume
V_0 = 0.4*w*h; % target volume

X = 0.4*ones(nel, 1); % initial density
X_filt = zeros(nel,1);
X_max = 1*ones(nel, 1);
X_min = 1e-3*ones(nel, 1);
N = findNeighbors(Ex,Ey,RADIUS);
%% Plot geometry before optimization
figure
axis equal
eldraw2(Ex, Ey, [1, 2, 0]);
hold on
plot(mean(Ex(42, :)), mean(Ey(42, :)), 'r*');
drawCircle(mean(Ex(42, :)), mean(Ey(42, :)), RADIUS, 1, 1);
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
tol = 1e-2;
X_norm = 1;
alpha = ALPHA;
%frames(loadSteps) = struct('cdata',[],'colormap',[]);
while X_norm > tol && c < 600
    for ie=1:nel
        X_filt(ie) = sum(N(ie,:)*X);
    end
    p = P_PAR;
    c = c + 1;
    K = sysK(k_0, X_filt.^p, Edof, ndof);
    K = sparse(K);
    u = solveq(K, F, bc);
    g_0plot = [g_0plot g_0(u)];
    X_next = optimization(X, u, Edof, k_0, V, V_0, X_max, X_min, p, alpha,X_filt,N);

    X_norm = norm(X_next - X);
    X = X_next;
    disp(c);
end
    for ie=1:nel
        X_filt(ie) = sum(N(ie,:)*X);
    end
%% Post-processing plotting
figure;
superdraw2(Ex, Ey, X);
axis equal;
annotation('textbox',...
    [0.14 0.94 0.34 0.06], ...
    'String', sprintf("Initial: %.1f; p = %d; Term. tol: %.1e; \\alpha = %.1f; Iterations: %d", ...
    0.4, p, tol, alpha, c), ...
    'FitBoxToText','on');

figure;
plot(g_0plot);
title("Objective function value");
ylabel("g_0");
xlabel("Iterations");
annotation('textbox',...
    [0.45 0.73 0.34 0.14], ...
    'String', sprintf("Initial: %.1f; p = %d; Term. tol: %.1e;\n \\alpha = %.1f; Iterations: %d", ...
    0.4, p, tol, alpha, c), ...
    'FitBoxToText','on');