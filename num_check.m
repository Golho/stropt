clearvars; close all;
%% OPTIONS
MESH_SCALE = 10;
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

X_start = VOLUME_PROCENT*ones(nel, 1); % initial density
X = X_start;
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
X_old1 = X;
X_old2 = X;
low = 0; 
upp = 1;
iter = 0; % counter
tol = 1e-2;
X_norm = 1;
alpha = ALPHA;
frames = []; %struct('cdata',[],'colormap',[]);
p = P_PAR;
%%
for ie=1:nel
    X_filt(ie) = sum(N(ie,:)*X);
end
K = sysK(k_0, X_filt.^p, Edof, ndof);
K = sparse(K);
if CASE == 2
    K(I ~= 0) = K(I ~= 0) + k;
end
u = solveq(K, F, bc);
g_00 = g_2(u);
%%
g_0diff = zeros(nel, 1);
h = 1e-8;
for e = 1:nel
    
    X(e) = X(e) + h;
%     for ie=1:nel
%             X_filt(ie) = sum(N(ie,:)*X);
%     end
    X_filt = X;
    iter = iter + 1;
    X_old2 = X_old1;
    X_old1 = X;

    K = sysK(k_0, X_filt.^p, Edof, ndof);
    K = sparse(K);
    if CASE == 2
        K(I ~= 0) = K(I ~= 0) + k;
    end
    u = solveq(K, F, bc);

    g_0diff(e) = g_2(u);
    X(e) = X_start(e);
end

gradienten = (g_0diff - g_00) / h;