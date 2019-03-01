clear all 
clc
E = 3e9;
v = 0.4;
V_max = 0.4*0.10*0.3;
[Edof,Enod,coord,Dof,F,I,bc]=genMesh(0.3,0.1,30,10);
nel = length(Enod(:,1));
X = zeros(nel,1);
X_old = zeros(nel,1);
K = zeros(length(F));
[ex,ey] = coordxtr(Edof,coord,Dof,4);
eldraw2(ex,ey,[2,4,1]);

a = zeros(length(F));
aold = zeros(length(F));
tol = 0.001;
rsd = 1;
D = [E/(1-v^2),v*E/(1-v^2),0
    v*E/(1-v^2),E/(1-v^2),0
    0,0,E*(1-v)/(2*(1-v^2))];

while rsd > tol
    K = zeros(length(F));
for ie= 1: nel
    elemdof = Edof(ie,2:end);
    Ke = planre([ex(ie,1),ex(ie,3)],[ey(ie,1),ey(ie,3)],[1,0.01],D);
    K(elemdof,elemdof) = K(elemdof,elemdof)+ Ke;
end
a = solveq(K,F,bc);
rsd = norm(aold-a);
aold=a;
end
eldisp2(ex,ey,a(Edof(:,2:end)),[1,4,1],1e4);
