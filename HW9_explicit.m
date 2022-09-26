% Solving of 2D Temperature eqn
% with finite differences


% 0 Clearing variables and figures
clear all % clearing memory
clf % clearing figures

% 1 Define numerical model
xsize = 100000; % horizontal, in m
ysize = 100000; % vertical, in m
nx=35; % grid resolution in horizontal direction
ny=45;% grid resolution in vertical direction
k=3; % diffusivity
dx=xsize/(nx-1); % horizontal grid step in m
dy=ysize/(ny-1); % vertical grid step in m
xp=-dx/2:dx:xsize+dx/2;
yp=-dy/2:dy:ysize+dy/2;

RHOCP=sparse(ny+1, nx+1);
kRHOCP=sparse(ny+1, nx+1); % this is a matrix containing k/RHOCP
T0=sparse(ny+1, nx+1);
for j=1:1:nx+1
    for i=1:1:ny+1
        radius = sqrt(((xsize+dx)/2-(j-1)*dx)^2 + ((ysize+dy)/2-dy*(i-1))^2);
        % sqrt(((xsize+dx)/2-xp(j))^2 + ((ysize+dy)/2-yp(i))^2); % distance from centre
        if (radius>= 20000)
            RHOCP(i, j) = 3300*1100; % density outside the radius of the plume
            T0(i, j)=1573;
        else
            RHOCP(i, j) = 3200*1000; % density inside the plume
            T0(i, j)=1873;
        end
        kRHOCP(i, j)=k/RHOCP(i, j);
    end
end

dt=dy^2/(4*max(max(kRHOCP))); % ny>nx => dy<dx => min is dy
% WHY NOT calculate 4*k/min(min(RHOCP))? 
% or min(dx,dy)^2*min(min(RHOCP))/(4k)

T=sparse(ny+1, nx+1);
for timestep=1:1:10
    for j=1:1:nx+1
        for i=1:1:ny+1

            if(i>2 && j>2 && i<ny && j<nx)
                T(i, j)=T0(i, j)+k*dt*((T0(i,j-1)-2*T0(i,j)+T0(i,j+1))/dx^2 + (T0(i-1,j)-2*T0(i,j)+T0(i+1,j))/dy^2)/RHOCP(i, j);
            end
            if(i<=2||j==2||i>=ny||j==nx)
             T(i, j)=1573;   
            end
            if (j==1)
                T(i, j)= T(i, j+1);
            end
            if (j==nx+1)
                T(i, j)= T(i, j-1); % if j==ny+1
            end
        end
    end
    T0=T;
end
% timesteps? BC?

disp(T(17, 15));

figure(1); clf
pcolor(xp, yp, T) ; colormap('Jet')
title('temperature T')
shading flat
colorbar
axis ij








