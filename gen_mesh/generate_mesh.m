%
%Matlab code to generate first mesh - DistMesh is used
%
addpath './distmesh'

n  =  11;
x_b = linspace(-0.4,0.4,n);
y_b = 0.0625*exp(-25*x_b.^2); %bump equation
vert = zeros(n+5,2); %total contour points
vert(2:n+1,1) = x_b;
vert(2:n+1,2) = y_b;
vert(1,1) = -1.5; % x right-top vertex
vert(1,2) = 0; % y right-top vertex
vert(n+2,1) = 1.5; % x right-top vertex
vert(n+2,2) = 0; % y right-top vertex
vert(n+3,1) = 1.5; % x right-top vertex
vert(n+3,2) = .8; % y right-top vertex
vert(n+4,1) = -1.5; % x left-top vertex
vert(n+4,2) = .8; % y left-top vertex
vert(n+5,1) = -1.5; % x left-bottom vertex
vert(n+5,2) = 0; % y left-bottom vertex
vert2 = zeros(n+5,2); %total contour points
vert2(2:n+1,1) = x_b;
vert2(2:n+1,2) = y_b;
vert2(1,1) = -1.5; % x right-top vertex
vert2(1,2) = 0; % y right-top vertex
vert2(n+2,1) = 1.5; % x right-top vertex
vert2(n+2,2) = 0; % y right-top vertex
vert2(n+3,1) = 1.5; % x right-top vertex
vert2(n+3,2) = -0.5; % y right-top vertex
vert2(n+4,1) = -1.5; % x left-top vertex
vert2(n+4,2) = -0.5; % y left-top vertex
vert2(n+5,1) = -1.5; % x left-bottom vertex
vert2(n+5,2) = 0; % y left-bottom ver
vert3 = zeros(n,2); %total contour points
vert3(:,1) = x_b;
vert3(:,2) = y_b;
fd =@(p) ddiff( drectangle(p,-1.5,1.5,0.,0.8), dpoly(p,vert2) );
fh=@(p) 0.05 + 0.5*dpoly(p,vert3);
[p,t]=distmesh2d(fd,fh,0.05,[-1.5,0; 1.5,.8],vert); %genereate mesh
size(t) %check number of elements

dlmwrite('nodes0',p, 'delimiter','\t') % store nodes coordinates
dlmwrite('element0',t, 'delimiter','\t') % store triangles points