% clear memory and figures
clear all % memory
clf % figures

% 1 define numerical model
xsize = 100000; % horizontal, m
ysize = 100000; % vertical, m
Nx=35; % grid resolution in horizontal direction
Ny=45;% grid resolution in vertical direction
dx=xsize/(Nx-1); % horizontal grid step in m
dy=ysize/(Ny-1); % vertical grid step in m
x=0:dx:xsize; % horizontal coordinates of grid points
y=0:dy:ysize; % vertical coordinates of grid points
gy=10; % defining gravitational constant
r=20000; % defining radius of plume
radius=sparse(Ny, Nx);
ETA=1e+19;

% defining the density depending on the location of the point
RHO = sparse(Ny, Nx);
for j=1:1:Nx
    for i=1:1:Ny
        radius (i, j) = sqrt((xsize/2-dx*(j-1))^2 + (ysize/2-dy*(i-1))^2); % distance from centre
        if (radius (i, j)> r)
            RHO(i, j) = 3300; % density outside the radius of the plume
        else RHO(i, j) = 3200; % density inside the plume
        end
    end
end


figure(1)
subplot (2, 1, 1); colormap('Jet')
pcolor(x, y, RHO) % color map

% 2 define global matrices L() and R()
N=Nx*Ny; % total number of unknowns
L=sparse(N, N); % coeff in the left side of equations
R=zeros(N,1); % values for the right hand side of eqns


% 3 composing global matrices
% going through all grid points
for j=1:1:Nx % from 1 with step 1 (goes 1 by 1 towards Nx)
    for i=1:1:Ny
        % define global index of the current equation
        kg=(j-1)*Ny+i;
        % decide on the type of current equation (BC or PE)
        if (j==1 || i==1 || j==Nx || i==Ny)
            % BC equation 1*omega(i,j)=0
            L(kg, kg)=1; % coeff of omega on the edges is 1 and the RHS is 0, there is just one global index and one local idex kg for each eqn/omega
            R(kg, 1)=0; % RHS
        else
            % Poisson eqn d2omega/dx2+ d2omega/dy2 = gy*(RHO(i, j+1)-RHO(i, j-1))/(2*ETA*dx)
            %        omega2
            %         |
            % omega1----omega3----omega5
            %         |
            %        omega4
            % (omega1-2omega3+omega5)/dx2 + (omega2-2omega3+omega4)/dy2 = gy*(RHO(i, j+1)-RHO(i, j-1))/(2*ETA*dx)
            % left hand side
            L(kg,kg-Ny)=1/dx^2; %omega1
            L(kg,kg-1)=1/dy^2; %omega2
            L(kg,kg)=-2/dx^2-2/dy^2; %omega3
            L(kg,kg+1)=1/dy^2; %omega4
            L(kg,kg+Ny)=1/dx^2; %omega5
            % RHS
            R(kg, 1)=gy*(RHO(i, j+1)-RHO(i, j-1))/(2*ETA*dx);
        end
    end
end

% 4 compute solution
S=L\R;

% 5 reload S into omega
% going through all grid points
omega=zeros(Ny, Nx); % define empty geometrical array
for j=1:1:Nx % from 1 with step 1 (goes 1 by 1 towards Nx)
    for i=1:1:Ny
        % define global index of the current equation S
        kg=(j-1)*Ny+i;
        % reload S(kg) => omega(i, j)
        omega(i, j)= S(kg);
    end
end



% visualisation of omega
figure(2)

subplot (2, 1, 1); colormap('Jet')
pcolor(x, y, omega) % color map
shading interp % smoothing colours
colorbar % give colorbar with numbers
title('omega(i, j) as colour map')

subplot (2, 1, 2); colormap('Jet')
surf (x, y, omega)
shading interp
light % shed light on the surface
lighting phong % make lighting nicer
colorbar
title('omega(i, j) as 3D surface')

% printing omega(7, 5) for reference values
disp(omega(7, 5));

% 6 composing matrices for XI
% define global matrices P() and Q()
P=sparse(N, N); % coeff in the left side of equations
Q=zeros(N,1); % values for the right hand side of eqns
% going through all grid points
for j=1:1:Nx % from 1 with step 1 (goes 1 by 1 towards Nx)
    for i=1:1:Ny
        % define global index of the current equation
        kg=(j-1)*Ny+i;
        if (j==1 || i==1 || j==Nx || i==Ny)
            % BC equation 1*omega(i,j)=0
            P(kg, kg)=1;
            Q(kg, 1)=0;
        else
            % LHS
            P(kg,kg-Ny)=1/dx^2; %XI1
            P(kg,kg-1)=1/dy^2; %XI2
            P(kg,kg)=-2/dx^2-2/dy^2; %XI3
            P(kg,kg+1)=1/dy^2; %XI4
            P(kg,kg+Ny)=1/dx^2; %XI5
            % RHS
            Q(kg, 1)=omega(i, j);
        end
    end
end


% 7 solving Sk=L\R
Sk=P\Q;


% 8 reloading Sk into XI
XI=zeros(Ny, Nx); % define empty geometrical array
for j=1:1:Nx % from 1 with step 1 (goes 1 by 1 towards Nx)
    for i=1:1:Ny
        % define global index of the current equation Sk
        kg=(j-1)*Ny+i;
        % reload Sk(kg) => XI(i, j)
        XI(i, j)= Sk(kg);
    end
end

% printing XI(7, 5) for reference values
disp(XI(7, 5));

% visualisation of XI
figure(3)

subplot (2, 1, 1); colormap('Jet')
pcolor(x, y, XI) % color map
shading interp % smoothing colours
colorbar % give colorbar with numbers
title('XI(i, j) as colour map')

subplot (2, 1, 2); colormap('Jet')
surf (x, y, XI)
shading interp
light % shed light on the surface
lighting phong % make lighting nicer
colorbar
title('XI(i, j) as 3D surface')


% 9 computing vx and vy for internal points

vx=zeros(Ny, Nx);
vy=zeros(Ny, Nx);
for j=1:1:Nx % from 1 with step 1 (goes 1 by 1 towards Nx)
    for i=1:1:Ny
        if (j==1 || i==1 || j==Nx || i==Ny)
            vx(i, j)=0;
            vy(i, j)=0;
        else
        vx(i, j)= (XI(i+1, j) - XI(i-1, j))/(2*dy);
        vy(i, j)= (XI(i, j-1) - XI(i, j+1))/(2*dx);
        end
    end
end


% visualisation of vx, vy
figure(4)
subplot (2, 1, 1); colormap('Jet')
pcolor(x, y, vx) % color map
title('vx(i, j)')
figure(5)
subplot (2, 1, 1); colormap('Jet')
pcolor(x, y, vy) % color map
title('vy(i, j)')

% printing vx(7, 5) for reference values
disp(vx(7, 5));

% printing vy(7, 5) for reference values
disp(vy(7, 5));
