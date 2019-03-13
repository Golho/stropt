clear all; clc; close all;
%% OPTIONS
PLOTNOW = 1;
FILTER = 1;
CASE = 2;
MOVIE = 0;

TERM_TOL = 1e-3; % termination tolerance
MESH_SCALE = 10;
P_PAR = 3; % p-parameter
ALPHA = 1.5; % alpha-parameter
RADIUS = 0.015; % Radius of Weight function
U_MAX = 5.0e-3; % Maximum displacement in case 2
VOLUME_PROCENT = 0.4;
MMA_S = [0.5, 1.2, 0.7];
DELTA = 1e-6; % Lower boundary of the density
%% Initializing parameters and generating mesh
% Material parameters
E = 3e9;
v = 0.4;
t = 0.01;
k = 1e7;
ptype = 1; % Plane stress
D = [E/(1-v^2),     v*E/(1-v^2),    0
     v*E/(1-v^2),    E/(1-v^2),      0
     0,              0,              E*(1-v)/(2*(1-v^2))];

% Generate mesh
w = 0.3;
h = 0.1;
nx = MESH_SCALE*3;
ny = MESH_SCALE*1;
[Edof, Enod, coord, Dof, F, I, bc] = genMesh(w, h, nx, ny, CASE);
P = 50;
F = sparse(F);
[Ex, Ey] = coordxtr(Edof,coord,Dof,4);
nbrElems = length(Enod(:,1));
nbrDofs = numel(Dof);

V_box = w*h*t;
V_max = VOLUME_PROCENT*V_box; % target volume
V_e = V_box / (nx*ny); % element volume

X = VOLUME_PROCENT*ones(nbrElems, 1); % initial density
X_min = zeros(nbrElems, 1);
X_max = 1*ones(nbrElems, 1);

% Create the weight matrix N
if FILTER
    N = findNeighbors(Ex,Ey,RADIUS);
else
    N = sparse(eye(nbrElems));
end
%% Plot geometry before optimization
figure
axis equal
eldraw2(Ex, Ey, [1, 1, 0]);
title("Mesh of the original geometry");
xlabel("x [m]");
ylabel("y [m]");
hold on
plot(mean(Ex(42, :)), mean(Ey(42, :)), 'r*');
drawCircle(mean(Ex(42, :)), mean(Ey(42, :)), RADIUS, 1, 1);

figGraph = figure;
subplot(2, 2, 1)
g0Plot = plot(0, 'LineWidth', 1);
title("g_0");
subplot(2, 2, 2);
g1Plot = plot(0, 'LineWidth', 1);
title("g_1");
subplot(2, 2, 3);
g2Plot = plot(0, 'LineWidth', 1);
title("g_2");

figStruct = figure;
if CASE == 1
    structPlot(1) = superdraw2(Ex, Ey, X);
    structPlot(2) = superdraw2(-Ex, Ey, X);
elseif CASE == 2
    structPlot(1) = superdraw2(Ex, Ey, X);
    structPlot(2) = superdraw2(Ex, -Ey, X);
end
axis equal;
caxis([0 1]);
%% objective function
if CASE == 1
    g_0 = @(u) F'*u;
elseif CASE == 2
    g_0 = @(u) -I'*u / U_MAX;
end
g_1 = @(X) V_e/V_max*sum(N*X) - 1;
g_2 = @(u) F'*u /(P * U_MAX) - 1;
gplot = zeros(3, 600);
%% k_0
k_0 = cell(nbrElems, 1);
for e = 1:nbrElems
    ex = Ex(e, :);
    ey = Ey(e, :);
    k_0{e} = planre(ex([1, 3]), ey([1, 3]), [ptype t], D);
end
%% Optimization preparation
X_old1 = X;
X_old2 = X;
low = 0; 
upp = 1;
iter = 0; % iteration counter

X_norm = 1;
alpha = ALPHA;
frames = []; %struct('cdata',[],'colormap',[]);
p = P_PAR;
%% Optimization loop
while X_norm > TERM_TOL && iter < 600
    iter = iter + 1;
    
    X_filt = N*X;
    E = (1-DELTA)*X_filt.^p + DELTA;
    K = sysK(k_0, E, Edof, nbrDofs);
    K = sparse(K);
    if CASE == 2
        K(I ~= 0) = K(I ~= 0) + k;
    end
    u = solveq(K, F, bc);
    gplot(:, iter) = [g_0(u); g_1(X); g_2(u)];
    
    dEdrho = (1-DELTA)*p*X_filt.^(p-1);
    g_1gradient = V_e/V_max * ones(nbrElems, 1);
    if CASE == 1
        b_k = filterGradient(dEdrho, u, u, k_0, Edof, N) ...
            .* X .^(alpha+1);
        X_star_fun = @(lambda) (b_k ./ (lambda*g_1gradient) ).^(1/(alpha+1));
        g_1OC = @(X_var) g_1(X) + g_1gradient'*(X_var - X);
        
        X_next = optimization(X_min, X_max, g_1OC, X_star_fun);
    elseif CASE == 2
        z = solveq(K, I, bc);

        g_0gradient = filterGradient(dEdrho, z, u, k_0, Edof, N) / U_MAX;
        g_2gradient = filterGradient(dEdrho, u, u, k_0, Edof, N) / (-P*U_MAX);

        [X_next, low, upp] = mma_solver(iter, X, X_old1, X_old2, g_0gradient, ...
            [g_1(X); g_2(u)], [g_1gradient'; g_2gradient'], low, upp, [nbrElems, ...
            MMA_S, DELTA]);
    end
    
    X_norm = norm(X_next - X);
    
    X_old2 = X_old1;
    X_old1 = X;
    X = X_next;
    
    disp(iter);
    disp(g_0(u));
    %% Intermediate plotting
    if PLOTNOW
        for e = 1:length(structPlot)
            structPlot(e).CData = 1 - N*X;
        end
        g0Plot.YData = gplot(1, 1:iter);
        g1Plot.YData = gplot(2, 1:iter);
        g2Plot.YData = gplot(3, 1:iter);
        drawnow
    end
    
    % Plot movie frames
    if MOVIE
        frames = [frames getframe(figStruct)];
    end
end
%% Post-processing plotting
figure(figStruct);
if CASE == 1
    superdraw2(Ex, Ey, N*X);
    superdraw2(-Ex, Ey, N*X);
elseif CASE == 2
    superdraw2(Ex, Ey, N*X);
    superdraw2(Ex, -Ey, N*X);
end
axis equal;
annotation('textbox',...
    [0.14 0.94 0.34 0.06], ...
    'String', sprintf("Initial: %.1f; p = %d; Term. tol: %.1e; \\alpha = %.1f; Iterations: %d", ...
    0.4, p, TERM_TOL, alpha, iter), ...
    'FitBoxToText','on');

figure;
plot(gplot(1, 1:iter-1));
title("Objective function value");
ylabel("g_0");
xlabel("Iterations");
annotation('textbox',...
    [0.45 0.73 0.34 0.14], ...
    'String', sprintf("Initial: %.1f; p = %d; Term. tol: %.1e;\n \\alpha = %.1f; Iterations: %d", ...
    0.4, p, TERM_TOL, alpha, iter), ...
    'FitBoxToText','on');

figure;
plot(gplot(2, 1:iter));
title("Volume constraint value");
ylabel("g_1");
xlabel("Iterations");
annotation('textbox',...
    [0.45 0.73 0.34 0.14], ...
    'String', sprintf("Initial: %.1f; p = %d; Term. tol: %.1e;\n \\alpha = %.1f; Iterations: %d", ...
    0.4, p, TERM_TOL, alpha, iter), ...
    'FitBoxToText','on');

figure;
plot(gplot(3, 1:iter));
title("Maximum displacement constraint value");
ylabel("g_2");
xlabel("Iterations");
annotation('textbox',...
    [0.45 0.73 0.34 0.14], ...
    'String', sprintf("Initial: %.1f; p = %d; Term. tol: %.1e;\n \\alpha = %.1f; Iterations: %d", ...
    0.4, p, TERM_TOL, alpha, iter), ...
    'FitBoxToText','on');