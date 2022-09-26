% solving of 2D Poisson equation using Gauss Seidel iterative method
% d2FI/dx2+ d2FI/dy2 = 1
% with finite differences

% clear memory and figures
clear all % memory
clf % figures

% 1 define numerical model
xsize = 1; % horizontal, m
ysize = 1; % vertical, m
Nx=35; % grid resolution in horizontal direction
Ny=45;% grid resolution in vertical direction
dx=xsize/(Nx-1); % horizontal grid step in m
dy=ysize/(Ny-1); % vertical grid step in m
x=0:dx:xsize; % horizontal coordinates of grid points
y=0:dy:ysize; % vertical coordinates of grid points

% 2 define initial approximation for FI(i, j)
FI = zeros(Ny, Nx);
FInew = sparse(Ny, Nx);


% 3 define relaxation parameter theta
theta = 1.5;



% 4 iteratice cycle
FInew = FI; % initial value for FInew
    for niter=1:1:10 % number of iterations
     for j=2:1:Nx-1
        for i=2:1:Ny-1 % did not include outer points
                       % going through all points of the grid
           
              FInew(i, j) = FI(i, j) + theta * ( 1 -((FI(i,j-1)-2*FI(i,j)+FI(i,j+1))/dx^2 + (FI(i-1,j)-2*FI(i,j)+FI(i+1,j))/dy^2))/(-2/dx^2-2/dy^2);
              FI(i, j)=FInew(i, j); % update of each unknown is done immediately after each iteration
        end
     end
    end % ending the loops


% 5 visualisation
figure(1)

subplot (2, 1, 1); colormap('Jet')
pcolor(x, y, FInew) % color map
shading interp % smoothing colours
colorbar % give colorbar with numbers
title('FInew(i, j) as colour map')

subplot (2, 1, 2); colormap('Jet')
surf (x, y, FInew)
shading interp
light % shed light on the surface
lighting phong % make lighting nicer
colorbar
title('FInew(i, j) as 3D surface')

disp(FInew(7, 5)); % printed reference value