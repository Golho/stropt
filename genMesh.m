function [Edof,Enod,coord,Dof,F,I,bc]=genMesh(w,h,nx,ny)
% PURPOSE
%  Generate FE-mesh for a rectangular design domain.
%    NOTE 1:  the user must manually modify the code, see Step 3.
%    NOTE 2:  this function is only able to handle 2D four-node quad
%             elements of equal side lengths. Mind the input arguments.
%
% INPUT:   h: height of the box
%          l: width of the box
%          nx: number of elements in the horizontal direction
%          ny: number of elements in the vertical direction
%
% OUTPUT:
%          Edof:  element DOF connectivity matrix. dim(Edof) = nelm x 9
%                 The first row of Edof contain element numbers.
%          Enod:  element nodal connectivity matrix. dim(Enod) = nelm x 5
%          coord: global coordinate matrix
%          Dof:   global DOF matrix. dim(Dof) = nnods x 2
%          F:     load vector. dim(F) = ndof x 1. F is nonzero in the DOF
%                 where the load is applied.
%          I:     "identification vector" whose entries are zero in all
%                 DOFs except in the DOF connected to the spring.
%

dofnod=2;  % Number of dofs in each node
nelm=nx*ny;% Number of elements

% --- Step 1: Establish coord
if( abs(w/nx-h/ny) > 1e-8)
    error('Not quadratic elements')
end
le=w/nx;
nodesx=nx+1; nodesy=ny+1;
nnods=nodesx*nodesy;
coord=zeros(nnods,2);
ycoord=0;xcoord=0; % start at origin
nod=0;
for i=1:nodesy
    for j=1:nodesx
        nod=nod+1;
        coord(nod,:)=[xcoord ycoord];
        xcoord=xcoord+le;
    end
    ycoord=ycoord+le;
    if(xcoord>=le)xcoord=0;end
end
% --- end step 1

% --- step 2: Enod, Dof and Edof
Enod=zeros(nelm,5);
nnx=nx+1;
el=0;
for i=1:ny
    for j=1:nx
        el=el+1;
        a=nnx*(i-1)+j; % Only valid for the node numbering: down left -> down right and then
        b=nnx*i+j;     % successiveley going "up" for each row
        Enod(el,:)=[el a a+1 b+1 b];
    end
end
% Dof matrix
nnods=(nx+1)*(ny+1);
Dof=zeros(nnods,dofnod);
for i=1:nnods
    Dof(i,:)=[i*2-1 i*2];
end
% Edof matrix
Edof=zeros(nelm,dofnod*4+1);
for ie=1:nelm
    for j=1:4
        nod=Enod(ie,j+1);
        xdof=nod*2-1;
        ydof=nod*2;
        Edof(ie,j*2)=xdof;
        Edof(ie,j*2+1)=ydof;
    end
    Edof(ie)=ie;
end
% --- end step 2

% extract element coordinates
[ex,ey]=coordxtr(Edof,coord,Dof,4);

% --- Step 3: boundary conditions creation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note that this part of the code is where the user
% needs to modify the code according to
% their needs and the specific problem to be solved.
% (comment/uncomment code below)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% -- Step 3.1 - symmetry at left side of the mesh (for the MBB beam)
% Find elements at left side for symmetry BC
% Check for 0 in ex

symbcEl=zeros(ny,1); % note: ny
sze=size(ex(1,:),2);
yie=1;
for ie=1:nelm
    for i=1:sze
        if(ex(ie,i)==0)
            symbcEl(yie)=ie;
            yie=yie+1;
            break;
        end
    end
end
% Now create the actual BC. Start with the symmetry conditions
%
bcsymdofs=zeros(ny+1,1); % should be equal to the number of ynods
for i=1:ny
    xsymdof1=Edof(symbcEl(i),2); % pick out the two xdofs from the given
    xsymdof2=Edof(symbcEl(i),8); % elements using Edof
    bcsymdofs(i)  = xsymdof1;
    bcsymdofs(i+1)= xsymdof2;
    % In this way we overwrite values for the bcsymdofs. But it is OK,
    % it is supposed to be that way to avoid copies of the ysymdofs.
    ydofF=xsymdof2+1; % The last y-dof is where the load should be applied
end
% Fixed y-dof at right bottom corner.
% (nx+1) is the node number there.
fixedbcdofs=[(nx+1)*2-1;(nx+1)*2];  

zv=zeros(length(bcsymdofs),1); % zero vector
bcsym=[bcsymdofs zv];          % symmetry bc
bcfix=[fixedbcdofs(2) 0];      % fixed bc 

bc=[bcsym;bcfix];
F=zeros(nnods*dofnod,1);       % external force vector
F(ydofF)=-1;
I = F*0.; % not used for this problem
% -- end Step 3.1
%
% -- Step 3.2 - symmetry at bottom side of the mesh (for the mechanism)
% Find elements at bottom side for symmetry BC
% Check for 0 in ey

% symbcEl=zeros(nx,1); % note: nx
% sze=size(ey(1,:),2);
% xie=1;
% for ie=1:nelm
%     for i=1:sze
%         if(ey(ie,i)==0)
%             symbcEl(xie)=ie;
%             xie=xie+1;
%             break;
%         end
%     end
% end
% % Now create the actual BC. Start with the symmetry conditions
% bcsymdofs=zeros(nx+1,1); % should be equal to the number of xnods
% for i=1:nx
%     ysymdof1=Edof(symbcEl(i),3); % pick out the two ydofs from the given elements
%     ysymdof2=Edof(symbcEl(i),5); % using Edof
%     bcsymdofs(i)  = ysymdof1;
%     bcsymdofs(i+1)= ysymdof2;
%     % In this way we overwrite values for the bcsymdofs. But it is OK,
%     % it is supposed to be that way to avoid copies of the ysymdofs.
% end
% %
% factor = 0.25; % percentage of the (upper part of the) left boundary
%                % that is clamped/fixed.
% fixedbcdofs=[];
% nn = max(max(Enod));
% tmp = zeros(nn,2);
% for ie=1:nelm
%     for i=1:4
%         % left boundary and [factor] top of it
%         if(ex(ie,i)==0 && (ey(ie,i) > (1.-factor)*max(max(ey))) )
%             nod=Enod(ie,i+1);
%             xdof=nod*2-1;
%             ydof=nod*2;
%             
%             tmp(nod,1)=nod;
%             if(tmp(nod,2) == 0) % no double
%                 tmp(nod,2) = 1;
%                 fixedbcdofs=[fixedbcdofs; xdof];
%                 fixedbcdofs=[fixedbcdofs; ydof];
%             end
%             
%         end
%     end
% end
% 
% zv=zeros(length(bcsymdofs),1); % zero vector
% bcsym=[bcsymdofs zv];          % symmetry bc
% bcfix=[fixedbcdofs'; zeros(length(fixedbcdofs),1)'];  % fixed bc
% bc=[bcfix bcsym']';
%   
% F=zeros(nnods*dofnod,1);       % external force vector
% F(max(max(Edof)))=-1;          % apply in top right y-DOF
% 
% I = zeros(nnods*dofnod,1);   % "identification vector":
% I(1) = -1;                   % nonzero in the one of the DOF 
%                              % connected to the spring element 
% % -- end Part 3.2

% END
