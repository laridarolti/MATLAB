% 0 clear memory and figures
clear all % memory
clf % figures


% 1.1 define numerical model
xsize = 100000; % horizontal, m
ysize = 100000; % vertical, m
Nx=35; % grid resolution in horizontal direction
Ny=45;% grid resolution in vertical direction
dx=xsize/(Nx-1); % horizontal grid step in m
dy=ysize/(Ny-1); % vertical grid step in m
x=0:dx:xsize+dx; % horizontal coordinates of grid points
y=0:dy:ysize+dy; % vertical coordinates of grid points
xp=-dx/2:dx:xsize+dx/2;
yp=-dy/2:dy:ysize+dy/2;
xvx=0:dx:xsize+dx;
yvx=-dy/2:dy:ysize+dy/2;
xvy=-dx/2:dx:xsize+dx/2;
yvy=0:dy:ysize+dy;
gy=10; % defining gravitational constant
r=20000; % defining radius of plume
radius=sparse(Ny, Nx);
ETA=1e+19;



% 1.2 defining the density depending on the location of the point
% density defined in vy points
RHO = sparse(Ny+1, Nx+1);
for j=1:1:Nx+1
    for i=1:1:Ny+1
        radius (i, j) = sqrt(((xsize+dx)/2-(j-1)*dx)^2 + ((ysize)/2-dy*(i-1))^2); % distance from centre
        if (radius (i, j)>= r)
            RHO(i, j) = 3300; % density outside the radius of the plume
        else RHO(i, j) = 3200; % density inside the plume
        end
    end
end



% 2 defining global matrices L and R

N=(Nx+1)*(Ny+1)*3; % total number of unknowns
L=sparse(N, N); % coeff in the left side of equations
R=zeros(N,1); % values for the right hand side of eqns


% 3 composing global matrices
for j=1:1:Nx+1
    for i=1:1:Ny+1 
        gp=((j-1)*(Ny+1)+i-1)*3+1;
        gvx=gp+1;
        gvy=gp+2;
        if (i==1 || j==1 || j>=Nx || i==Ny+1)
            L(gvx, gvx)=1;
            R(gvx, 1)=0;  
        else
            L(gvx, gvx-3*(Ny+1))=ETA/dx^2;
            L(gvx, gvx-3)=ETA/dy^2;
            L(gvx, gvx)=-2*ETA*(1/dx^2+1/dy^2);
            L(gvx, gvx+3)=ETA/dy^2;
            L(gvx, gvx+3*(Ny+1))=ETA/dx^2;
            L(gvx, gp)=1/dx;
            L(gvx, gp+3*(Ny+1))=-1/dx;
            R(gvx, 1)=0;
        end
        if (i==1 || j==1 || j==Nx+1 || i>=Ny)
            L(gvy, gvy)=1;
            R(gvy, 1)=0;
        else
            L(gvy, gvy-3*(Ny+1))=ETA/dx^2;
            L(gvy, gvy-3)=ETA/dy^2;
            L(gvy, gvy)=-2*ETA*(1/dx^2+1/dy^2);
            L(gvy, gvy+3)=ETA/dy^2;
            L(gvy, gvy+3*(Ny+1))=ETA/dx^2;
            L(gvy, gp)=1/dy;
            L(gvy, gp+3)=-1/dy;
            R(gvy, 1)=-RHO(i, j)*gy;
        end
        
        if (i==1 || j==1 || j==Nx+1 || i==Ny+1 ||(i==2 && j==2))
            L(gp, gp)=1;
            R(gp, 1)=0;

        else
            L(gp, gvx-3*(Ny+1))=-1/dx;
            L(gp, gvx)=1/dx;
            L(gp, gvy-3)=-1/dy;
            L(gp, gvy)=1/dy;
            R(gp, 1)=0;
        end
        end
end


% 4 solve matrices
S=L\R;


% 5 reload solution into vx, vy, p
p=zeros(Ny+1, Nx+1);
vx=zeros(Ny+1, Nx+1);
vy=zeros(Ny+1, Nx+1);  % define empty geometrical arrays
for j=1:1:Nx+1 % from 1 with step 1 
    for i=1:1:Ny+1
        % define global index of the current equation S
        gp=((j-1)*(Ny+1)+i-1)*3+1;
        gvx=gp+1;
        gvy=gp+2;
        % reload S(kg) => p, vx, vy(i, j)
        p(i, j)= S(gp);
        vx(i, j)= S(gvx);
        vy(i, j)= S(gvy);
    end
end


figure(1); clf
pcolor(xvy, yvy, RHO);
colormap('Jet')
shading flat
colorbar

figure(2); clf
pcolor(xp, yp, p) ; colormap('Jet')
shading flat
colorbar
axis ij

figure(3); clf
pcolor(xvx, yvx, vx) ; colormap('Jet')
shading flat
colorbar
axis ij

figure(4); clf
pcolor(xvy, yvy, vy) ; colormap('Jet')
shading flat
colorbar
axis ij

% reference values
disp(p(7, 5));
disp(vx(7, 5));
disp(vy(7, 5));

