% 0 clear memory and figures
clear all % memory
clf % figures


% DEFINE NUMERICAL MODEL
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
BC_type=-1; % free slip=-1, no slip =1
RHOmelt=2500;
ETA_melt=10;



% RHO, K_PHI_VY defined in vy points
RHO = zeros(Ny+1, Nx+1);
K_phi_vy= zeros(Ny+1, Nx+1);
K_phi_vx= zeros(Ny+1, Nx+1);
ETA_phi= zeros(Ny+1, Nx+1);
ETA_B= zeros(Ny+1, Nx+1);
ETA_P= zeros(Ny+1, Nx+1);

for j=1:1:Nx+1
    for i=1:1:Ny+1
        radius1= sqrt(((xsize+dx)/2-(j-1)*dx)^2 + ((ysize)/2-dy*(i-1))^2); % radius for vy pts
        radius2= sqrt(((xsize+dx)/2-(j-1)*dx)^2 + ((ysize+dy)/2-dy*(i-1))^2); % radius for p pts
        radius3= sqrt(((xsize)/2-(j-1)*dx)^2 + ((ysize+dy)/2-dy*(i-1))^2);  % radius for vx pts
        radius4= sqrt((xsize/2-(j-1)*dx)^2 + (ysize/2-dy*(i-1))^2);  % radius for  basic nodal pts
        if (radius1>= r)
            RHO(i, j) = 3300; % density mantle
            K_phi_vy(i, j)=1e-12;
        else
            RHO(i, j) = 3200; % density inside the plume
            K_phi_vy(i, j)=1e-11;
        end
        if (radius2>= r)
            ETA_P(i, j) = 1E19; % viscosity outside the radius of the plume
            ETA_phi(i, j)=1e21;
        else
            ETA_P(i, j) = 1E18; % viscosity inside the plume
            ETA_phi(i, j)=1e20;
        end
        if (radius3>= r)
            K_phi_vx(i, j)=1e-12;
        else
            K_phi_vx(i, j)=1e-11;
        end
        if (radius4>= r)
            ETA_B(i, j) = 1E19; % viscosity outside the radius of the plume
        else
            ETA_B(i, j) = 1E18; % viscosity inside the plume
            %disp(ETA_B(i, j));
        end
    end
end


% 2 defining global matrices L and R

N=(Nx+1)*(Ny+1)*6; % total number of unknowns
L=sparse(N, N); % coeff in the left side of equations
R=zeros(N,1); % values for the right hand side of eqns


% 3 COMPOSING GLOBAL MATRICES
for j=1:1:Nx+1
    for i=1:1:Ny+1
        gp=((j-1)*(Ny+1)+i-1)*6+1;
        gvx=gp+1;
        gvy=gp+2;
        gpmelt=gp+3;
        gvxD=gp+4;
        gvyD=gp+5;

        % CONTINUITY EQUATION SOLID
        if (i==1 || j==1 || j==Nx+1 || i==Ny+1)
            L(gp, gp)=1;
            R(gp, 1)=0;
        elseif (i==2 && j==2)
            L(gp, gp)=1;
            R(gp, 1)=3.3E9;
        else
            L(gp, gvx-6*(Ny+1))=-1/dx;
            L(gp, gvx)=1/dx;
            L(gp, gvy-6)=-1/dy;
            L(gp, gvy)=1/dy;

            % NEW PART (p-pmelt)/eta_phi
            L(gp,gp)=1/ETA_phi(i, j);
            L(gp, gpmelt)=-1/ETA_phi(i, j);
            R(gp, 1)=0;
        end

        % x-stokes equation - bulk
        if (j==1 || j>=Nx)
            L(gvx, gvx)=1;
            R(gvx, 1)=0;
        elseif (i==1)
            L(gvx, gvx)=1;
            L(gvx, gvx+6)=BC_type;
            R(gvx, 1)=0;
        elseif (i==Ny+1)
            L(gvx, gvx)=1;
            L(gvx, gvx-6)=BC_type;
            R(gvx, 1)=0;
        else
            % vx velocities
            L(gvx, gvx-6*(Ny+1))=2*ETA_P(i, j)/dx^2;
            L(gvx, gvx-6)=ETA_B(i-1, j)/dy^2;
            L(gvx, gvx)=-2*(ETA_P(i, j)+ETA_P(i, j+1))/dx^2-(ETA_B(i-1, j)+ETA_B(i, j))/dy^2;
            L(gvx, gvx+6)=ETA_B(i, j)/dy^2;
            L(gvx, gvx+6*(Ny+1))=2*ETA_P(i, j+1)/dx^2;
            % vy velocities
            L(gvx, gvy-6)=ETA_B(i-1, j)/(dx*dy);
            L(gvx, gvy)=-ETA_B(i, j)/(dx*dy);
            L(gvx, gvy+6*Ny)=-ETA_B(i-1, j)/(dx*dy);
            L(gvx, gvy+6*(Ny+1))=ETA_B(i, j)/(dx*dy);
            % pressure
            L(gvx, gp)=1/dx;
            L(gvx, gp+6*(Ny+1))=-1/dx;
            % RHS
            R(gvx, 1)=0;
        end

        % y-stokes equation - bulk
        if (i==1 || i>=Ny)
            L(gvy, gvy)=1;
            R(gvy, 1)=0;
        elseif (j==1)
            L(gvy, gvy)=1;
            L(gvy, gvy+6*(Ny+1))=BC_type;
            R(gvy, 1)=0;
        elseif (j==Nx+1)
            L(gvy, gvy)=1;
            L(gvy, gvy-6*(Ny+1))=BC_type;
            R(gvy, 1)=0;
        else
            % vx velocities
            L(gvy, gvx-6*(Ny+1))=ETA_B(i, j-1)/(dx*dy);
            L(gvy, gvx-6*Ny)=-ETA_B(i, j-1)/(dx*dy);
            L(gvy, gvx)=-ETA_B(i, j)/(dx*dy);
            L(gvy, gvx+6)=ETA_B(i, j)/(dx*dy);
            % vy velocities
            L(gvy, gvy-6*(Ny+1))=ETA_B(i, j-1)/dx^2;
            L(gvy, gvy-6)=2*ETA_P(i, j)/dy^2;
            L(gvy, gvy)=-2*(ETA_P(i, j)+ETA_P(i+1, j))/dy^2-(ETA_B(i, j-1)+ETA_B(i, j))/dx^2;
            L(gvy, gvy+6)=2*ETA_P(i+1, j)/dy^2;
            L(gvy, gvy+6*(Ny+1))=ETA_B(i, j)/dx^2;
            % pressure
            L(gvy, gp)=1/dy;
            L(gvy, gp+6)=-1/dy;
            %RHS
            R(gvy, 1)=-RHO(i, j)*gy;
        end

        % CONTINUITY EQUATION MELT
        if(i==1 || j==1 || j==Nx+1 || i==Ny+1)
            L(gpmelt, gpmelt)=1;
            R(gpmelt, 1)=0;
        else
            L(gpmelt, gvxD-6*(Ny+1))=-1/dx;
            L(gpmelt, gvxD)=1/dx;
            L(gpmelt, gvyD-6)=-1/dy;
            L(gpmelt, gvyD)=1/dy;

            % NEW PART (p-pmelt)/eta_phi
            L(gpmelt,gp)=-1/ETA_phi(i, j);
            L(gpmelt, gpmelt)=1/ETA_phi(i, j);
            R(gpmelt, 1)=0;
        end

        % x-Darsy equation
        if (j==1 || j>=Nx)
            L(gvxD, gvxD)=1;
            R(gvxD, 1)=0;
        elseif (i==1)
            L(gvxD, gvxD)=1;
            L(gvxD, gvxD+6)=BC_type;
            R(gvxD, 1)=0;
        elseif (i==Ny+1)
            L(gvxD, gvxD)=1;
            L(gvxD, gvxD-6)=BC_type;
            R(gvxD, 1)=0;
        else
            % vxD velocities
            L(gvxD, gvxD)=ETA_melt/K_phi_vx(i, j);
            % pressure
            L(gvxD, gpmelt)=-1/dx;
            L(gvxD, gpmelt+6*(Ny+1))=1/dx;
            % RHS
            R(gvxD, 1)=0;
        end

        % y-Darsy equation - melt
        if (i==1 || i>=Ny)
            L(gvyD, gvyD)=1;
            R(gvyD, 1)=0;
        elseif (j==1)
            L(gvyD, gvyD)=1;
            L(gvyD, gvyD+6*(Ny+1))=BC_type;
            R(gvyD, 1)=0;
        elseif (j==Nx+1)
            L(gvyD, gvyD)=1;
            L(gvyD, gvyD-6*(Ny+1))=BC_type;
            R(gvyD, 1)=0;
        else
            % vyD velocities
            L(gvyD, gvyD)=ETA_melt/K_phi_vy(i, j);
            % pressure
            L(gvyD, gpmelt)=-1/dy;
            L(gvyD, gpmelt+6)=1/dy;
            %RHS
            R(gvyD, 1)=RHOmelt*gy;
        end
    end
end


% 4 solve matrices
S=L\R;


% 5 reload solution into vx, vy, p
p=zeros(Ny+1, Nx+1);
pmelt=zeros(Ny+1, Nx+1);
vx=zeros(Ny+1, Nx+1);
vxD=zeros(Ny+1, Nx+1);
vy=zeros(Ny+1, Nx+1);
vyD=zeros(Ny+1, Nx+1);% define empty geometrical arrays
for j=1:1:Nx+1 % from 1 with step 1
    for i=1:1:Ny+1
        % define global index of the current equation S
        gp=((j-1)*(Ny+1)+i-1)*6+1;
        gvx=gp+1;
        gvy=gp+2;
        gpmelt=gp+3;
        gvxD=gp+4;
        gvyD=gp+5;
        % reload S(gp) => p, vx, vy, pmelt, vxD, vyD(i, j)
        p(i, j)= S(gp);
        pmelt(i, j)=S(gpmelt);
        vx(i, j)= S(gvx);
        vxD(i, j)=S(gvxD);
        vy(i, j)= S(gvy);
        vyD(i, j)=S(gvyD);
    end
end


figure(1); clf
pcolor(xvy, yvy, RHO);
title('density RHO')
colormap('Jet')
shading flat
colorbar

figure(2); clf
pcolor(xvy, yvy, K_phi_vy);
title('K_phi_vy')
colormap('Jet')
shading flat
colorbar

figure(3); clf
pcolor(xvx, yvx, K_phi_vx);
title('K_phi_vx')
colormap('Jet')
shading flat
colorbar

figure(4); clf
pcolor(xp, yp, ETA_P) ; colormap('Jet')
title('ETA_P viscosity in pressure points')
shading flat
colorbar
axis ij

figure(5); clf
pcolor(x, y, ETA_B) ; colormap('Jet')
title('ETA_B viscosity in basic nodal points')
shading flat
colorbar
axis ij

figure(6); clf
pcolor(xp, yp, ETA_phi) ; colormap('Jet')
title('ETA_phi viscosity in basic nodal points')
shading flat
colorbar
axis ij


figure(7); clf
pcolor(xp, yp, p) ; colormap('Jet')
title('pressure p')
shading flat
colorbar
axis ij

figure(8); clf
pcolor(xvx, yvx, vx) ; colormap('Jet')
title('vx velocity')
shading flat
colorbar
axis ij

figure(9); clf
pcolor(xvy, yvy, vy) ; colormap('Jet')
title('vy velocity')
shading flat
colorbar
axis ij

figure(10); clf
pcolor(xp, yp, pmelt) ; colormap('Jet')
title('pressure pmelt')
shading flat
colorbar
axis ij

figure(11); clf
pcolor(xvx, yvx, vxD) ; colormap('Jet')
title('vxD velocity')
shading flat
colorbar
axis ij

figure(12); clf
pcolor(xvy, yvy, vyD) ; colormap('Jet')
title('vyD velocity')
shading flat
colorbar
axis ij

figure(13); clf
pcolor(xp, yp,p- pmelt) ; colormap('Jet')
title('pressure difference p- pmelt')
shading flat
colorbar
axis ij

% reference values
disp(p(7, 5));
disp(vx(7, 5));
disp(vy(7, 5));
disp(pmelt(7, 5));
disp(vxD(7, 5));
disp(vyD(7, 5));