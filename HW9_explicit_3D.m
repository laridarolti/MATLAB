% Solving of 2D Temperature eqn
% with finite differences


% 0 Clearing variables and figures
clear all % clearing memory
clf % clearing figures

% 1 Define numerical model
xsize = 100000; % horizontal, in m
ysize = 100000; % vertical, in m
zsize = 100000;
nx=35; % grid resolution in horizontal direction
ny=45;% grid resolution in vertical direction
nz=25;
k=3; % diffusivity
dx=xsize/(nx-1); % horizontal grid step in m
dy=ysize/(ny-1); % vertical grid step in m
dz=zsize/(nz-1);
xp=-dx/2:dx:xsize+dx/2;
yp=-dy/2:dy:ysize+dy/2;
zp=-dz/2:dz:zsize+dz/2;

RHOCP=zeros(ny+1, nx+1, nz+1);
kRHOCP=zeros(ny+1, nx+1, nz+1); % this is a matrix containing k/RHOCP
T0=zeros(ny+1, nx+1, nz+1);
for l=1:1:nz+1
    for i=1:1:ny+1
        for j=1:1:nx+1
            radius = sqrt(((xsize+dx)/2-(j-1)*dx)^2 + ((ysize+dy)/2-dy*(i-1))^2+((zsize+dz)/2-dz*(l-1))^2);
            if (radius>= 20000)
                RHOCP(i, j, l) = 3300*1100; % density outside the radius of the plume
                T0(i, j, l)=1573;
            else
                RHOCP(i, j, l) = 3200*1000; % density inside the plume
                T0(i, j, l)=1873;
            end
            kRHOCP(i, j, l)=k/RHOCP(i, j, l);
        end
    end
end

dt=dy^2/(4*max(max(max(kRHOCP)))); % ny>nx => dy<dx => min is dy
% WHY NOT calculate 4*k/min(min(RHOCP))? 
% or min(dx,dy)^2*min(min(RHOCP))/(4k)

T=zeros(ny+1, nx+1, nz+1);
for timestep=1:1:10
    for l=1:1:nz+1
        for j=1:1:nx+1
            for i=1:1:ny+1

                if(i>2 && j>2 && i<ny && j<nx && l>2 && l<nz)
                    T(i, j, l)=T0(i, j, l)+k*dt*((T0(i,j-1, l)-2*T0(i,j, l)+T0(i,j+1, l))/dx^2 + (T0(i-1,j, l)-2*T0(i,j, l)+T0(i+1,j, l))/dy^2+ (T0(i,j, l-1)-2*T0(i,j, l)+T0(i,j, l+1))/dz^2)/RHOCP(i, j, l);
                end
                if(i<=2||j==2||i>=ny||j==nx || l<=2 || l>=nz)
                    T(i, j, l)=1573;
                end
                if (j==1)
                    T(i, j, l)= T(i, j+1,l);
                end
                if (j==nx+1)
                    T(i, j, l)= T(i, j-1, l); % if j==ny+1
                end
            end
        end
    end
    T0=T;
end
% timesteps? BC?

disp(T(17, 15, 13));

% figure(1); clf
% pcolor(xp, yp, zp, T) ; colormap('Jet')
% title('temperature T')
% shading flat
% colorbar
% axis ij








