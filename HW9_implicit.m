% 0 clear memory and figures
clear all % memory
clf % figures


% 1.1 define numerical model
xsize = 100000; % horizontal, in m
ysize = 100000; % vertical, in m
nx=35; % grid resolution in horizontal direction
ny=45;% grid resolution in vertical direction
k=3; % diffusivity
dx=xsize/(nx-1); % horizontal grid step in m
dy=ysize/(ny-1); % vertical grid step in m
xp=-dx/2:dx:xsize+dx/2;
yp=-dy/2:dy:ysize+dy/2;

% 2 Define sparse matrix L() and vector R()
N=(nx+1)*(ny+1)*3; % total number of unknowns
L=sparse(N, N); % coeff in the left side of equations
R=zeros(N,1); % values for the right hand side of eqns


RHOCP=sparse(ny+1, nx+1);
kRHOCP=sparse(ny+1, nx+1); % this is a matrix containing k/RHOCP
T=sparse(ny+1, nx+1);
for j=1:1:nx+1
    for i=1:1:ny+1
        radius = sqrt(((xsize+dx)/2-(j-1)*dx)^2 + ((ysize+dy)/2-dy*(i-1))^2);
        % sqrt(((xsize+dx)/2-xp(j))^2 + ((ysize+dy)/2-yp(i))^2); % distance from centre
        if (radius>= 20000)
            RHOCP(i, j) = 3300*1100; % density outside the radius of the plume
            T(i, j)=1573;
        else
            RHOCP(i, j) = 3200*1000; % density inside the plume
            T(i, j)=1873;
        end
        kRHOCP(i, j)=k/RHOCP(i, j);
    end
end

dt=dy^2/(4*max(max(kRHOCP))); % ny>nx => dy<dx => min is dy



for timestep=1:1:10


    % 3 composing global matrices
    for j=1:1:nx+1
        for i=1:1:ny+1
            gp=((j-1)*(ny+1)+i-1)*3+1;


            % temperature equation
            if(i<=2||j==2||i>=ny||j==nx)
                L(gp, gp)=1;
                R(gp, 1)=1573;
            elseif (j==1)
                L(gp, gp)=1;
                L(gp, gp+3*(ny+1))=-1;
                R(gp, 1)=0;
            elseif (j==nx+1)
                L(gp, gp)=1;
                L(gp, gp-3*(ny+1))=-1;
                R(gp, 1)=0;
            else
                L(gp, gp-3*(ny+1))=-dt*k/dx^2;
                L(gp, gp)=RHOCP(i, j)+2*dt*k/dy^2+2*dt*k/dx^2;
                L(gp, gp-3)=-dt*k/dy^2;
                L(gp, gp+3*(ny+1))=-dt*k/dx^2;
                L(gp, gp+3)=-dt*k/dy^2;
                R(gp, 1)=RHOCP(i, j)*T(i, j);
            end

        end
    end


    % 4 solve matrices
    S=L\R;


    % 5 reload solution into vx, vy, p
    T=zeros(ny+1, nx+1);
    for j=1:1:nx+1 % from 1 with step 1
        for i=1:1:ny+1
            % define global index of the current equation S
            gp=((j-1)*(ny+1)+i-1)*3+1;
            % reload S(kg) => p, vx, vy(i, j)
            T(i, j)= S(gp);
        end
    end
end

figure(1); clf
pcolor(xp, yp, T) ; colormap('Jet')
title('temperature T')
shading flat
colorbar
axis ij

disp(T(17, 15));