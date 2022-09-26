% 0 clear memory and figures
clear all % memory
clf % figures


% 1.1 define numerical model
xsize = 100000; % horizontal, m
ysize = 100000; % vertical, m
nx=35; % grid resolution in horizontal direction
ny=45;% grid resolution in vertical direction
dx=xsize/(nx-1); % horizontal grid step in m
dy=ysize/(ny-1); % vertical grid step in m
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
BC_type=-1; %  free slip



% define marker arrays
nxm=(nx-1)*5; % number of markers in horizontal direction
nym=(ny-1)*5;
nm=nxm*nym; % total number of markers
xm=zeros(nm, 1); % x coordinates of markers
ym=zeros(nm, 1);
RHOm=zeros(nm, 1); % density value for markers in kg/m3
ETAm=zeros(nm, 1); % viscosity values for markers






% defining initial positions of markers and value of RHOm and ETAm for them

dxm=xsize/nxm; % average horizontal distance between markers in m
dym=ysize/nym;
k=1;
for jm=1:1:nxm
    for im=1:1:nym
        % coordinates
        xm(k)=dxm/2+(jm-1)*dxm ; %+(rand-0.5)*dxm; %horizontal
        ym(k)=dym/2+(im-1)*dym ; %+(rand-0.5)*dym; %vertical
% 
%         if (xm(k)>x(nx))
%             xm(k)=xm(k)-xsize;
%         end
%         if (xm(k)<x(1))
%             xm(k)=xm(k)+xsize;
%         end
%         if (ym(k)>y(ny))
%             ym(k)=ym(k)-ysize;
%         end
%         if (ym(k)<y(1))
%             ym(k)=ym(k)+ysize;
%         end
        % update marker counter
        k=k+1;
    end
end

% 2 defining some matrices

N=(nx+1)*(ny+1)*3; % total number of unknowns
L=sparse(N, N); % coeff in the left side of equations
R=zeros(N,1); % values for the right hand side of eqns
RHO = sparse(ny+1, nx+1);
ETA_P = sparse(ny+1, nx+1);
ETA_B = sparse(ny+1, nx+1);



% time stepping
for timestep=1:1:10
    %define RHOm(m) and ETAm(m)
    for m=1:1:nm
        radius = sqrt((xm(m)-(xsize)/2)^2 + (ym(m)-(ysize)/2)^2); % distance from centre
        if (radius>=r)
            RHOm(m) = 3300; % density outside the radius of the plume
            ETAm(m) = 1E19; % viscosity outside the plume
        else
            RHOm(m) = 3200; % density inside the plume
            ETAm(m) = 1E18;
        end
    end

    % initial values for sum arrays for ETAP, ETAB, RHO
    ETA_P_SUM=zeros(ny+1, nx+1);
    ETA_P_WT=zeros(ny+1, nx+1);
    ETA_B_SUM=zeros(ny+1, nx+1);
    ETA_B_WT=zeros(ny+1, nx+1);
    RHO_SUM=zeros(ny+1, nx+1);
    RHO_WT=zeros(ny+1, nx+1);


    % calculate ETA_P

    for m=1:1:nm

        % indices
        i=fix(((ym(m)-yp(1))/dy))+1;
        j=fix(((xm(m)-xp(1))/dx))+1;


        % conditions for i, j
        if (j<1)
            j=1;
        end
        if (i<1)
            i=1;
        end
        if(j>nx-1)
            j=nx-1;
        end
        if(i>ny-1)
            i=ny-1;
        end

        % weights for pressure pts
        dxmj=(xm(m)-xp(j));
        dymi=(ym(m)-yp(i));

        wtmij=(1-dxmj/dx)*(1-dymi/dy);
        wtmi1j=(1-dxmj/dx)*dymi/dy;
        wtmij1=dxmj/dx*(1-dymi/dy);
        wtmi1j1=dxmj/dx*dymi/dy;

        % calculating numerator and denominator for ETA_P
        ETA_P_SUM(i, j)=ETA_P_SUM(i, j)+wtmij*ETAm(m);
        ETA_P_WT(i, j)=ETA_P_WT(i, j)+wtmij;
        ETA_P_SUM(i+1, j)=ETA_P_SUM(i+1, j)+wtmi1j*ETAm(m);
        ETA_P_WT(i+1, j)=ETA_P_WT(i+1, j)+wtmi1j;
        ETA_P_SUM(i, j+1)=ETA_P_SUM(i, j+1)+wtmij1*ETAm(m);
        ETA_P_WT(i, j+1)=ETA_P_WT(i, j+1)+wtmij1;
        ETA_P_SUM(i+1, j+1)=ETA_P_SUM(i+1, j+1)+wtmi1j1*ETAm(m);
        ETA_P_WT(i+1, j+1)=ETA_P_WT(i+1, j+1)+wtmi1j1;

        % indices
        i= fix(((ym(m)-y(1)))/dy)+1;
        j= fix(((xm(m)-x(1)))/dx)+1;

        % conditions for i, j
        if (j<1)
            j=1;
        end
        if (i<1)
            i=1;
        end
        if(j>nx-1)
            j=nx-1;
        end
        if(i>ny-1)
            i=ny-1;
        end

        % weights for basic nodal pts, use same variables, previous values
        % not needed
        dxmj=(xm(m)-x(j));
        dymi=(ym(m)-y(i));

        wtmij=(1-dxmj/dx)*(1-dymi/dy);
        wtmi1j=(1-dxmj/dx)*(dymi/dy);
        wtmij1=dxmj/dx*(1-dymi/dy);
        wtmi1j1=dxmj/dx*(dymi/dy);

        % calculating numerator and denominator for ETA_B
        ETA_B_SUM(i, j)=ETA_B_SUM(i, j)+wtmij*ETAm(m);
        ETA_B_WT(i, j)=ETA_B_WT(i, j)+wtmij;
        ETA_B_SUM(i+1, j)=ETA_B_SUM(i+1, j)+wtmi1j*ETAm(m);
        ETA_B_WT(i+1, j)=ETA_B_WT(i+1, j)+wtmi1j;
        ETA_B_SUM(i, j+1)=ETA_B_SUM(i, j+1)+wtmij1*ETAm(m);
        ETA_B_WT(i, j+1)=ETA_B_WT(i, j+1)+wtmij1;
        ETA_B_SUM(i+1, j+1)=ETA_B_SUM(i+1, j+1)+wtmi1j1*ETAm(m);
        ETA_B_WT(i+1, j+1)=ETA_B_WT(i+1, j+1)+wtmi1j1;

        % indices
        i=fix(((ym(m)-yvy(1)))/dy)+1;
        j=fix(((xm(m)-xvy(1)))/dx)+1;

        % conditions for i, j
        if (j<1)
            j=1;
        end
        if (i<1)
            i=1;
        end
        if(j>nx-1)
            j=nx-1;
        end
        if(i>ny-1)
            i=ny-1;
        end

        % weights for vy pts, use same variables, previous values
        % not needed
        dxmj=abs(xm(m)-xvy(j));
        dymi=abs(ym(m)-yvy(i));
        
        wtmij=(1-dxmj/dx)*(1-dymi/dy);
        wtmi1j=(1-dxmj/dx)*dymi/dy;
        wtmij1=dxmj/dx*(1-dymi/dy);
        wtmi1j1=dxmj/dx*dymi/dy;

        % calculating numerator and denominator for RHO
        RHO_SUM(i, j)=RHO_SUM(i, j)+wtmij*RHOm(m);
        RHO_WT(i, j)=RHO_WT(i, j)+wtmij;
        RHO_SUM(i+1, j)=RHO_SUM(i+1, j)+wtmi1j*RHOm(m);
        RHO_WT(i+1, j)=RHO_WT(i+1, j)+wtmi1j;
        RHO_SUM(i, j+1)=RHO_SUM(i, j+1)+wtmij1*RHOm(m);
        RHO_WT(i, j+1)=RHO_WT(i, j+1)+wtmij1;
        RHO_SUM(i+1, j+1)=RHO_SUM(i+1, j+1)+wtmi1j1*RHOm(m);
        RHO_WT(i+1, j+1)=RHO_WT(i+1, j+1)+wtmi1j1;
    end

    % calculating ETAP, ETAB, RHO
    for j=1:1:nx+1
        for i=1:1:ny+1
            if(ETA_P_WT(i, j)>0)
                ETA_P(i, j)=ETA_P_SUM(i, j)/ETA_P_WT(i, j);
            else
                ETA_P(i, j)=ETA_P(i, j);
            end
            if (ETA_B_WT(i, j)>0)
                ETA_B(i, j)=ETA_B_SUM(i, j)/ETA_B_WT(i, j);
            else
                ETA_B(i, j)=ETA_B(i, j);
            end
            if (RHO_WT(i, j)>0)
                RHO(i, j)=RHO_SUM(i, j)/RHO_WT(i, j);
            else
                RHO(i, j)=RHO(i, j);
            end
        end
    end




    % 3 composing global matrices
    for j=1:1:nx+1
        for i=1:1:ny+1
            gp=((j-1)*(ny+1)+i-1)*3+1;
            gvx=gp+1;
            gvy=gp+2;


            % continuity equation
            if (i==1 || j==1 || j==nx+1 || i==ny+1)
                L(gp, gp)=1;
                R(gp, 1)=0;
            elseif (i==2 && j==2)
                L(gp, gp)=1;
                R(gp, 1)=3.3E9;
            else
                L(gp, gvx-3*(ny+1))=-1/dx;
                L(gp, gvx)=1/dx;
                L(gp, gvy-3)=-1/dy;
                L(gp, gvy)=1/dy;
                R(gp, 1)=0;
            end

            % x-stokes equation

            if (j==1 || j>=nx)
                L(gvx, gvx)=1;
                R(gvx, 1)=0;
            elseif (i==1)
                L(gvx, gvx)=1;
                L(gvx, gvx+3)=BC_type;
                R(gvx, 1)=0;
            elseif (i==ny+1)
                L(gvx, gvx)=1;
                L(gvx, gvx-3)=BC_type;
                R(gvx, 1)=0;
            else
                % vx velocities
                L(gvx, gvx-3*(ny+1))=2*ETA_P(i, j)/dx^2;
                L(gvx, gvx-3)=ETA_B(i-1, j)/dy^2;
                L(gvx, gvx)=-2*(ETA_P(i, j)+ETA_P(i, j+1))/dx^2-(ETA_B(i-1, j)+ETA_B(i, j))/dy^2;
                L(gvx, gvx+3)=ETA_B(i, j)/dy^2;
                L(gvx, gvx+3*(ny+1))=2*ETA_P(i, j+1)/dx^2;

                % vy velocities
                L(gvx, gvy-3)=ETA_B(i-1, j)/(dx*dy);
                L(gvx, gvy)=-ETA_B(i, j)/(dx*dy);
                L(gvx, gvy+3*ny)=-ETA_B(i-1, j)/(dx*dy);
                L(gvx, gvy+3*(ny+1))=ETA_B(i, j)/(dx*dy);

                % pressure
                L(gvx, gp)=1/dx;
                L(gvx, gp+3*(ny+1))=-1/dx;

                % RHS
                R(gvx, 1)=0;
            end

            % y-stokes equation
            if (i==1 || i>=ny)
                L(gvy, gvy)=1;
                R(gvy, 1)=0;
            elseif (j==1)
                L(gvy, gvy)=1;
                L(gvy, gvy+3*(ny+1))=BC_type;
                R(gvy, 1)=0;
            elseif (j==nx+1)
                L(gvy, gvy)=1;
                L(gvy, gvy-3*(ny+1))=BC_type;
                R(gvy, 1)=0;
            else

                % vx velocities
                L(gvy, gvx-3*(ny+1))=ETA_B(i, j-1)/(dx*dy);
                L(gvy, gvx-3*ny)=-ETA_B(i, j-1)/(dx*dy);
                L(gvy, gvx)=-ETA_B(i, j)/(dx*dy);
                L(gvy, gvx+3)=ETA_B(i, j)/(dx*dy);

                % vy velocities
                L(gvy, gvy-3*(ny+1))=ETA_B(i, j-1)/dx^2;
                L(gvy, gvy-3)=2*ETA_P(i, j)/dy^2;
                L(gvy, gvy)=-2*(ETA_P(i, j)+ETA_P(i+1, j))/dy^2-(ETA_B(i, j-1)+ETA_B(i, j))/dx^2;
                L(gvy, gvy+3)=2*ETA_P(i+1, j)/dy^2;
                L(gvy, gvy+3*(ny+1))=ETA_B(i, j)/dx^2;

                % pressure
                L(gvy, gp)=1/dy;
                L(gvy, gp+3)=-1/dy;

                %RHS
                R(gvy, 1)=-RHO(i, j)*gy;
            end
        end
    end


    % 4 solve matrices
    S=L\R;


    % 5 reload solution into vx, vy, p
    p=zeros(ny+1, nx+1);
    vx=zeros(ny+1, nx+1);
    vy=zeros(ny+1, nx+1);
    for j=1:1:nx+1 % from 1 with step 1
        for i=1:1:ny+1
            % define global index of the current equation S
            gp=((j-1)*(ny+1)+i-1)*3+1;
            gvx=gp+1;
            gvy=gp+2;
            % reload S(kg) => p, vx, vy(i, j)
            p(i, j)= S(gp);
            vx(i, j)= S(gvx);
            vy(i, j)= S(gvy);
        end
    end



    % calculating dt
    % defining dt
    dt=1e+30;
    vxmax=max(max(abs(vx)));
    vymax=max(max(abs(vy)));
    if (dt*vxmax>dx*0.1)
        dt=dx*0.1/vxmax;
    end
    if (dt*vymax>dy*0.1)
        dt=dy*0.1/vymax;
    end




    % interpolationg vx to vxm
    vxm=sparse(nm);
    vym=sparse(nm);
    for m=1:1:nm
        i=fix(((ym(m)-yvx(1)))/dy)+1;
        j=fix(((xm(m)-xvx(1)))/dx)+1;
        if (j<1)
            j=1;
        end
        if (i<1)
            i=1;
        end
        if(j>nx-1)
            j=nx-1;
        end
        if(i>ny-1)
            i=ny-1;
        end

        dxmj=(xm(m)-xvx(j));
        dymi=(ym(m)-yvx(i));

        wtmij=(1-dxmj/dx)*(1-dymi/dy);
        wtmi1j=(1-dxmj/dx)*dymi/dy;
        wtmij1=dxmj/dx*(1-dymi/dy);
        wtmi1j1=dxmj/dx*dymi/dy;


        vxm(m)=vx(i, j)*wtmij + vx(i+1, j)*wtmi1j+vx(i, j+1)*wtmij1+vx(i+1, j+1)*wtmi1j1;
        vym(m)=vy(i, j)*wtmij + vy(i+1, j)*wtmi1j+vy(i, j+1)*wtmij1+vy(i+1, j+1)*wtmi1j1;



        % moving markers
        xm(m)=xm(m)+vxm(m)*dt;
        ym(m)=ym(m)+vym(m)*dt;

        if (xm(m)>x(nx))
            xm(m)=xm(m)-xsize;
        end
        if (xm(m)<x(1))
            xm(m)=xm(m)+xsize;
        end
        if (ym(m)>y(ny))
            ym(m)=ym(m)-ysize;
        end
        if (ym(m)<y(1))
            ym(m)=ym(m)+ysize;
        end

    end

end

figure(1); clf
pcolor(xvy, yvy, RHO);
title('density RHO')
colormap('Jet')
shading flat
colorbar

figure(2); clf
pcolor(xp, yp, ETA_P) ; colormap('Jet')
title('ETA_P viscosity in pressure points')
shading flat
colorbar
axis ij

figure(3); clf
pcolor(x, y, ETA_B) ; colormap('Jet')
title('ETA_B viscosity in basic nodal points')
shading flat
colorbar
axis ij

figure(4); clf
pcolor(xp, yp, p) ; colormap('Jet')
title('pressure p')
shading flat
colorbar
axis ij

figure(5); clf
pcolor(xvx, yvx, vx) ; colormap('Jet')
title('vx velocity')
shading flat
colorbar
axis ij

figure(6); clf
pcolor(xvy, yvy, vy) ; colormap('Jet')
title('vy velocity')
shading flat
colorbar
axis ij

disp(p(7, 5));
disp(vx(7, 5));
disp(vy(7, 5));
