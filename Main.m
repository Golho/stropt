clear all; clc; close all;
%% OPTIONS
MESH_SCALE = 20;
P_PAR = 3; % p-parameter
ALPHA = 1.5; % alpha-parameter
RADIUS = 0.015; % Radius of Weight function
CASE = 2;
U_MAX = 0.5e-3; % Maximum displacement in case 2
VOLUME_PROCENT = 0.4;
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
nel = length(Enod(:,1));
ndof = numel(Dof);

V = w*h/(nx*ny); % element volume
V_0 = VOLUME_PROCENT*w*h; % target volume

X = VOLUME_PROCENT*ones(nel, 1); % initial density
X_filt = zeros(nel,1);
X_max = 1*ones(nel, 1);
X_min = 1e-2*ones(nel, 1);
N = findNeighbors(Ex,Ey,RADIUS);
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
%% objective function
if CASE == 1
    g_0 = @(u) F'*u;
elseif CASE == 2
    g_0 = @(u) -I'*u / U_MAX;
end
g_1 = @(X) V*sum(X)/V_0 - 1;
g_2 = @(u) F'*u /(P * U_MAX) - 1;
gplot = zeros(3, 600);
%% k_0
k_0 = cell(nel, 1);
for e = 1:nel
    ex = Ex(e, :);
    ey = Ey(e, :);
    k_0{e} = planre(ex([1, 3]), ey([1, 3]), [ptype t], D);
end
%% Optimization loop
X_old1 = X;%zeros(nel, 1);
X_old2 = X;%zeros(nel, 1);
low = 0; 
upp = 1;
iter = 0; % counter
tol = 1e-3;
X_norm = 1;
alpha = ALPHA;
frames = []; %struct('cdata',[],'colormap',[]);
p = P_PAR;
while X_norm > tol && iter < 600
    iter = iter + 1;
    for ie=1:nel
        X_filt(ie) = sum(N(ie,:)*X);
    end
    
    % OBS; INGEN FILTRERING SKER
    X_filt = X;

    K = sysK(k_0, X_filt.^p, Edof, ndof);
    K = sparse(K);
    if CASE == 2
        K(I ~= 0) = K(I ~= 0) + k;
    end
    u = solveq(K, F, bc);
    gplot(:, iter) = [g_0(u); g_1(X); g_2(u)];
    if CASE == 1
        X_next = optimization(X, u, Edof, k_0, V, V_0, X_max, X_min, p, alpha,X_filt,N);
    elseif CASE == 2
        z = solveq(K, I, bc);
        fGrad = filterGradient(X_filt, p, k_0, u, Edof, N, 2);
        g_0gradient = zeros(nel, 1);
        g_2gradient = zeros(nel, 1);
        for e = 1:nel
            idof = Edof(e, 2:end);
            dK = sparse(ndof, ndof, 64);
            dK(idof, idof) = p*X(e)^(p-1)*k_0{e};
            g_0gradient(e) = u'*dK*z / U_MAX;
            g_2gradient(e) = -u'*dK*u / (P*U_MAX);
        end
        
        %g_0gradient = fGrad' * z / U_MAX;
        g_1gradient = V/V_0*ones(nel, 1);
        %g_2gradient = u*%-filterGradient(X_filt, p, k_0, u, Edof, N, 1) / ( P * U_MAX );

        [X_next, low, upp] = mma_solver(iter, X, X_old1, X_old2, g_0gradient, ...
            [g_1(X); g_2(u)], [g_1gradient, g_2gradient], low, upp, [nel, ...
            0.5, 0.5, 1.5, X_min(1)]);
    end
    X_norm = norm(X_next - X);
    
    X_old2 = X_old1;
    X_old1 = X;
    
    X = X_next;
    
    
    disp(iter);
%     if mod(iter, 5) == 0
%         figStructure = figure('visible', 'off');
%         superdraw2(Ex, Ey, X);
%         axis equal;
%         frames = [frames getframe(figStructure)];
%         close(figStructure);
%     end
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
    0.4, p, tol, alpha, iter), ...
    'FitBoxToText','on');

figure;
plot(gplot(1, 1:iter));
title("Objective function value");
ylabel("g_0");
xlabel("Iterations");
annotation('textbox',...
    [0.45 0.73 0.34 0.14], ...
    'String', sprintf("Initial: %.1f; p = %d; Term. tol: %.1e;\n \\alpha = %.1f; Iterations: %d", ...
    0.4, p, tol, alpha, iter), ...
    'FitBoxToText','on');

figure;
plot(gplot(2, 1:iter));
title("Volume constraint value");
ylabel("g_1");
xlabel("Iterations");
annotation('textbox',...
    [0.45 0.73 0.34 0.14], ...
    'String', sprintf("Initial: %.1f; p = %d; Term. tol: %.1e;\n \\alpha = %.1f; Iterations: %d", ...
    0.4, p, tol, alpha, iter), ...
    'FitBoxToText','on');

figure;
plot(gplot(3, 1:iter));
title("Maximum displacement constraint value");
ylabel("g_2");
xlabel("Iterations");
annotation('textbox',...
    [0.45 0.73 0.34 0.14], ...
    'String', sprintf("Initial: %.1f; p = %d; Term. tol: %.1e;\n \\alpha = %.1f; Iterations: %d", ...
    0.4, p, tol, alpha, iter), ...
    'FitBoxToText','on');