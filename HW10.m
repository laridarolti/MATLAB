% 0 clear memory and figures
clear all % memory
clf % figures


% 1.1 define numerical model
xsize = 100000; % horizontal, in m
ysize = 100000; % vertical, in m
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
BC_type=-1; % free slip
BC_pressure = 1e+5;
gy=10;
dtthermal= min(dx,dy)^2/(4*max([4/3350/800 3/3300/1000  2/3200/1100]));
timesum=0;

% FOR TESTCASES
nsteps=1; dt=6e+10; dispmax=0.5;


% define marker arrays
nxm=(nx-1)*4; % number of markers in horizontal direction
nym=(ny-1)*4;
nm=nxm*nym; % total number of markers
xm=zeros(nm, 1); % x coordinates of markers
ym=zeros(nm, 1);
dxm=xsize/nxm; % average horizontal distance between markers in m
dym=ysize/nym;
RHOm=zeros(nm, 1); % density value for markers in kg/m3
ETAm=zeros(nm, 1); % viscosity values for markers
CPm=zeros(nm, 1);
RHOCPm=zeros(nm, 1);
km=zeros(nm, 1);
kappam=zeros(nm, 1);
Tm=zeros(nm, 1);

m=1;
for jm=1:1:nxm
    for im=1:1:nym
        % coordinates
        xm(m)=dxm/2+(jm-1)*dxm ; %+(rand-0.5)*dxm; %horizontal
        ym(m)=dym/2+(im-1)*dym ; %+(rand-0.5)*dym; %vertical

        % STICKY AIR
        if (ym(m)<0.2*ysize)
            RHOm(m) = 1; % sticky air
            CPm(m)=3E6;
            km(m)=3000;
            ETAm(m) = 1E17;
            Tm(m)=273;
        end

        % LITHOSPHERE
        if (ym(m)>=0.2*ysize && ym(m)<0.4*ysize )
            RHOm(m) = 3350;
            CPm(m)=800;
            km(m)=4;
            ETAm(m) = 1E22;
            Tm(m)=-1027+ym(m)*13/200;
        end

        if (ym(m)>=0.4*ysize)

            % PLUME
            if (sqrt((xm(m)-50000)^2+(ym(m)-70000)^2)<20000)
                RHOm(m) = 3200;
                CPm(m)=1100;
                km(m)=2;
                ETAm(m) = 1E18;
                Tm(m)=1873;

                % MANTLE
            else
                RHOm(m) = 3300;
                CPm(m)=1000;
                km(m)=3;
                ETAm(m) = 1E19;
                Tm(m)=1573;
            end
        end
        RHOCPm(m)=RHOm(m)*CPm(m);
        m=m+1;
    end
end

% 2 defining some matrices

N=(nx+1)*(ny+1)*3; % total number of unknowns
L=sparse(N, N); % coeff in the left side of equations
R=zeros(N,1); % values for the right hand side of eqns
RHO = zeros(ny+1, nx+1);
RHOCP = zeros(ny+1, nx+1);
ETA_P = zeros(ny+1, nx+1);
ETA_B = zeros(ny+1, nx+1);
kx= zeros(ny+1, nx+1);
ky= zeros(ny+1, nx+1);
T0= zeros(ny+1, nx+1);




for timestep=1:1:nsteps

    % initial values for sum arrays for ETAP, ETAB, RHO, kx, ky
    ETA_P_SUM=zeros(ny+1, nx+1);
    ETA_P_WT=zeros(ny+1, nx+1);
    ETA_B_SUM=zeros(ny+1, nx+1);
    ETA_B_WT=zeros(ny+1, nx+1);
    RHO_SUM=zeros(ny+1, nx+1);
    RHO_WT=zeros(ny+1, nx+1);
    RHOCP_SUM=zeros(ny+1, nx+1);
    RHOCP_WT=zeros(ny+1, nx+1);
    kx_SUM= zeros(ny+1, nx+1);
    kx_WT= zeros(ny+1, nx+1);
    ky_SUM= zeros(ny+1, nx+1);
    ky_WT= zeros(ny+1, nx+1);
    T_SUM= zeros(ny+1, nx+1);
    T_WT= zeros(ny+1, nx+1);


    % calculate ETA_P, ETA_B, RHO (or RHOvy in vy points)
    for m=1:1:nm

        % I. ETA_P, T and RHOCP
        % indices
        i=fix(((ym(m)-yp(1))/dy))+1;
        j=fix(((xm(m)-xp(1))/dx))+1;

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

        % calculating numerator and denominator for RHOCP
        RHOCP_SUM(i, j)=RHOCP_SUM(i, j)+wtmij*RHOCPm(m);
        RHOCP_WT(i, j)=RHOCP_WT(i, j)+wtmij;
        RHOCP_SUM(i+1, j)=RHOCP_SUM(i+1, j)+wtmi1j*RHOCPm(m);
        RHOCP_WT(i+1, j)=RHOCP_WT(i+1, j)+wtmi1j;
        RHOCP_SUM(i, j+1)=RHOCP_SUM(i, j+1)+wtmij1*RHOCPm(m);
        RHOCP_WT(i, j+1)=RHOCP_WT(i, j+1)+wtmij1;
        RHOCP_SUM(i+1, j+1)=RHOCP_SUM(i+1, j+1)+wtmi1j1*RHOCPm(m);
        RHOCP_WT(i+1, j+1)=RHOCP_WT(i+1, j+1)+wtmi1j1;

        % calculating numerator and denominator for T
        T_SUM(i, j)=T_SUM(i, j)+wtmij*Tm(m)*RHOCPm(m);
        T_WT(i, j)=T_WT(i, j)+wtmij*RHOCPm(m);
        T_SUM(i+1, j)=T_SUM(i+1, j)+wtmi1j*Tm(m)*RHOCPm(m);
        T_WT(i+1, j)=T_WT(i+1, j)+wtmi1j*RHOCPm(m);
        T_SUM(i, j+1)=T_SUM(i, j+1)+wtmij1*Tm(m)*RHOCPm(m);
        T_WT(i, j+1)=T_WT(i, j+1)+wtmij1*RHOCPm(m);
        T_SUM(i+1, j+1)=T_SUM(i+1, j+1)+wtmi1j1*Tm(m)*RHOCPm(m);
        T_WT(i+1, j+1)=T_WT(i+1, j+1)+wtmi1j1*RHOCPm(m);


        % II. ETA_B
        % indices
        i= fix(((ym(m)-y(1)))/dy)+1;
        j= fix(((xm(m)-x(1)))/dx)+1;

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

        % III. RHOvy, ky
        % indices
        i=fix(((ym(m)-yvy(1)))/dy)+1;
        j=fix(((xm(m)-xvy(1)))/dx)+1;

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

        % calculating numerator and denominator for ky
        ky_SUM(i, j)=ky_SUM(i, j)+wtmij*km(m);
        ky_WT(i, j)=ky_WT(i, j)+wtmij;
        ky_SUM(i+1, j)=ky_SUM(i+1, j)+wtmi1j*km(m);
        ky_WT(i+1, j)=ky_WT(i+1, j)+wtmi1j;
        ky_SUM(i, j+1)=ky_SUM(i, j+1)+wtmij1*km(m);
        ky_WT(i, j+1)=ky_WT(i, j+1)+wtmij1;
        ky_SUM(i+1, j+1)=ky_SUM(i+1, j+1)+wtmi1j1*km(m);
        ky_WT(i+1, j+1)=ky_WT(i+1, j+1)+wtmi1j1;



        % IV. kx
        % indices
        i=fix(((ym(m)-yvx(1)))/dy)+1;
        j=fix(((xm(m)-xvx(1)))/dx)+1;

        % weights for vx pts, use same variables, previous values
        % not needed
        dxmj=abs(xm(m)-xvx(j));
        dymi=abs(ym(m)-yvx(i));

        wtmij=(1-dxmj/dx)*(1-dymi/dy);
        wtmi1j=(1-dxmj/dx)*dymi/dy;
        wtmij1=dxmj/dx*(1-dymi/dy);
        wtmi1j1=dxmj/dx*dymi/dy;

        % calculating numerator and denominator for kx
        kx_SUM(i, j)=kx_SUM(i, j)+wtmij*km(m);
        kx_WT(i, j)=kx_WT(i, j)+wtmij;
        kx_SUM(i+1, j)=kx_SUM(i+1, j)+wtmi1j*km(m);
        kx_WT(i+1, j)=kx_WT(i+1, j)+wtmi1j;
        kx_SUM(i, j+1)=kx_SUM(i, j+1)+wtmij1*km(m);
        kx_WT(i, j+1)=kx_WT(i, j+1)+wtmij1;
        kx_SUM(i+1, j+1)=kx_SUM(i+1, j+1)+wtmi1j1*km(m);
        kx_WT(i+1, j+1)=kx_WT(i+1, j+1)+wtmi1j1;

    end

    % calculating ETAP, ETAB, RHO, kx, ky, T0, RHOCP
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
            if (RHOCP_WT(i, j)>0)
                RHOCP(i, j)=RHOCP_SUM(i, j)/RHOCP_WT(i, j);
            else
                RHOCP(i, j)=RHOCP(i, j);
            end
            if (kx_WT(i, j)>0)
                kx(i, j)=kx_SUM(i, j)/kx_WT(i, j);
            else
                kx(i, j)=kx(i, j);
            end
            if (ky_WT(i, j)>0)
                ky(i, j)=ky_SUM(i, j)/ky_WT(i, j);
            else
                ky(i, j)=ky(i, j);
            end
            if (T_WT(i, j)>0)
                T0(i, j)=T_SUM(i, j)/T_WT(i, j);
            else
                T0(i, j)=T0(i, j);
            end
        end
    end
    T0(2:ny,1)=T0(2:ny,2); % Left BC: dT/dx=0
    T0(2:ny,nx+1)=T0(2:ny,nx); % Right BC: dT/dx=0
    T0(1,:)=2*273-T0(2,:); % Top BC: T=273
    T0(ny+1,:)=2*1573-T0(ny,:); % Bottom BC: T=1573



    % 3 composing global matrices L, R
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
                R(gp, 1)=BC_pressure;
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
                L(gvy, gvx-3*(ny+1))=ETA_B(i, j-1)/(dx*dy)-0.25*gy*dt*(RHO(i, j+1)-RHO(i, j-1))/(2*dx);
                L(gvy, gvx-3*ny)=-ETA_B(i, j-1)/(dx*dy)-0.25*gy*dt*(RHO(i, j+1)-RHO(i, j-1))/(2*dx);
                L(gvy, gvx)=-ETA_B(i, j)/(dx*dy)-0.25*gy*dt*(RHO(i, j+1)-RHO(i, j-1))/(2*dx);
                L(gvy, gvx+3)=ETA_B(i, j)/(dx*dy)-0.25*gy*dt*(RHO(i, j+1)-RHO(i, j-1))/(2*dx);

                % vy velocities
                L(gvy, gvy-3*(ny+1))=ETA_B(i, j-1)/dx^2;
                L(gvy, gvy-3)=2*ETA_P(i, j)/dy^2;
                L(gvy, gvy)=-2*(ETA_P(i, j)+ETA_P(i+1, j))/dy^2-(ETA_B(i, j-1)+ETA_B(i, j))/dx^2-gy*dt*(RHO(i+1, j)-RHO(i-1, j))/(2*dy);
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
    if (dt*vxmax>dx*dispmax)
        dt=dx*dispmax/vxmax;
    end
    if (dt*vymax>dy*dispmax)
        dt=dy*dispmax/vymax;
    end


    % interpolationg vx to vxm
    vxm=zeros(1,nm);
    vym=zeros(1,nm);
    vxm_eff=zeros(1,nm);
    vym_eff=zeros(1,nm);



    for m=1:1:nm
        xa=xm(m); ya=ym(m);
        for RK=1:1:4
            coeff=1;% coefficient for vxA, B, C, D in RK scheme
            if (RK==2 || RK==3)
                coeff=2; % changing coeff for B and C
            end

            i=fix(((ya-yvx(1)))/dy)+1;
            j=fix(((xa-xvx(1)))/dx)+1;
            dxmj=(xa-xvx(j));
            dymi=(ya-yvx(i));


            % NEW PART - BILIN INTERP
            vx13=(1-dxmj/dx)*vx(i,   j)+dxmj/dx*vx(i,   j+1);
            vx24=(1-dxmj/dx)*vx(i+1, j)+dxmj/dx*vx(i+1, j+1);

            vx13_corr_r= 0;
            vx24_corr_r= 0;
            vx13_corr_l= 0;
            vx24_corr_l= 0;

            if(j<nx-1)
                vx13_corr_r=0.5*(dxmj/dx-0.5)^2*(vx(i,   j)-2*vx(i,   j+1)+vx(i,   j+2));
                vx24_corr_r=0.5*(dxmj/dx-0.5)^2*(vx(i+1, j)-2*vx(i+1, j+1)+vx(i+1, j+2));
            end
            if (j>=2)
                vx13_corr_l=0.5*(dxmj/dx-0.5)^2*(vx(i,   j-1)-2*vx(i,   j)+vx(i,   j+1));
                vx24_corr_l=0.5*(dxmj/dx-0.5)^2*(vx(i+1, j-1)-2*vx(i+1, j)+vx(i+1, j+1));
            end

            if(xa>xvx(j)+0.5*dx) % CORRECTION TO RIGHT, pt to the right of pressure pt
                vx13=vx13+vx13_corr_r;
                vx24=vx24+vx24_corr_r;
            end
            if(xa<xvx(j)+0.5*dx) % CORRECTION TO LEFT, pt to the left of pressure pt
                vx13=vx13+vx13_corr_l;
                vx24=vx24+vx24_corr_l;
            end

            vxa=vx13*(1-dymi/dy)+ vx24*dymi/dy; % this is basically vxm(m), but do not want to change vxm(m) so used different variable

            i=fix(((ya-yvy(1)))/dy)+1;
            j=fix(((xa-xvy(1)))/dx)+1;
            dxmj=(xa-xvy(j));
            dymi=(ya-yvy(i));


            % NEW PART - BILIN INTERP
            vy13=(1-dymi/dy)*vy(i, j  )+(dymi/dy)*vy(i+1, j);
            vy24=(1-dymi/dy)*vy(i, j+1)+(dymi/dy)*vy(i+1, j+1);

            vy13_corr_r= 0;
            vy24_corr_r= 0;
            vy13_corr_l= 0;
            vy24_corr_l= 0;

            if(i<ny-1)
                vy13_corr_r=0.5*(dymi/dy-0.5)^2*(vy(i, j  )-2*vy(i+1, j  )+vy(i+2, j  ));
                vy24_corr_r=0.5*(dymi/dy-0.5)^2*(vy(i, j+1)-2*vy(i+1, j+1)+vy(i+2, j+1));
            end
            if(i>=2)
                vy13_corr_l=0.5*(dymi/dy-0.5)^2*(vy(i-1, j  )-2*vy(i, j  )+vy(i+1, j  ));
                vy24_corr_l=0.5*(dymi/dy-0.5)^2*(vy(i-1, j+1)-2*vy(i, j+1)+vy(i+1, j+1));
            end

            if(ya>yvy(i)+0.5*dy)
                vy13=vy13+vy13_corr_r;
                vy24=vy24+vy24_corr_r;
            end
            if(ya<yvy(i)+0.5*dy)
                vy13=vy13+vy13_corr_l;
                vy24=vy24+vy24_corr_l;
            end
            vya=vy13*(1-dxmj/dx)+ vy24*dxmj/dx;

            if (RK==1)
                vxm(m)=vxa; vym(m)=vya; % vxm and vym are the same as vxa and vya when rk = 1, as if there was no RK scheme
            end
            % CALCULATING EFFECTIVE VELOCITIES, UPDATING FOR PTS A, B, C, D
            vxm_eff(m)=vxm_eff(m) + coeff*vxa; % UPDATING vxm

            coeff1=1;
            if (RK==1 || RK==2)
                coeff1=2;
            end
            xa=xm(m)+vxa*dt/coeff1;

            vym_eff(m)=vym_eff(m) + coeff*vya;
            ya=ym(m)+vya*dt/coeff1;

        end

        % DIVIDING BY 6 AT THE END
        vxm_eff(m)=vxm_eff(m)/6;
        vym_eff(m)=vym_eff(m)/6;

        % moving markers
        xm(m)=xm(m)+vxm_eff(m)*dt;
        ym(m)=ym(m)+vym_eff(m)*dt;


    end

    NT=(nx+1)*(ny+1)*1; % total number of unknowns
    LT=sparse(NT, NT); % coeff in the left side of equations
    RT=zeros(NT,1); % values for the right hand side of eqns

    dt=min(dt, dtthermal);


    % 3 composing global matrices
    for j=1:1:nx+1
        for i=1:1:ny+1
            gp=((j-1)*(ny+1)+i-1)*1+1;
            if(i==1)
                LT(gp, gp)=1;
                LT(gp, gp+1)=1;
                RT(gp, 1)=2*273; % temp should be 273 at the top
            elseif(i==ny+1)
                LT(gp, gp)=1;
                LT(gp, gp-1)=1;
                RT(gp, 1)=2*1573;
            elseif (j==1)
                LT(gp, gp)=1;
                LT(gp, gp+1*(ny+1))=-1;
                RT(gp, 1)=0;
            elseif (j==nx+1)
                LT(gp, gp)=1;
                LT(gp, gp-1*(ny+1))=-1;
                RT(gp, 1)=0;
            else
                LT(gp, gp-1*(ny+1))=-dt*kx(i, j-1)/dx^2;
                LT(gp, gp-1)=-dt*ky(i-1, j)/dy^2;
                LT(gp, gp)=RHOCP(i, j)+dt*(ky(i, j)+ky(i-1, j))/dy^2+dt*(kx(i, j-1)+kx(i, j))/dx^2;
                LT(gp, gp+1)=-dt*ky(i, j)/dy^2;
                LT(gp, gp+1*(ny+1))=-dt*kx(i, j)/dx^2;
                RT(gp, 1)=RHOCP(i, j)*T0(i, j);
            end
        end
    end


    % 4 solve matrices
    ST=LT\RT;


    % 5 reload solution into T
    T=zeros(ny+1, nx+1);
    for j=1:1:nx+1 % from 1 with step 1
        for i=1:1:ny+1
            gp=((j-1)*(ny+1)+i-1)*1+1;
            % reload ST(gp) => T(i, j)
            T(i, j)= ST(gp);
        end
    end
    DT=T-T0;
    % MOVING T MARKERS
    if (nsteps==1)
        for m=1:1:nm
            i=fix(((ym(m)-yp(1)))/dy)+1;
            j=fix(((xm(m)-xp(1)))/dx)+1;

            dxmj=(xm(m)-xp(j));
            dymi=(ym(m)-yp(i));

            wtmij=(1-dxmj/dx)*(1-dymi/dy);
            wtmi1j=(1-dxmj/dx)*dymi/dy;
            wtmij1=dxmj/dx*(1-dymi/dy);
            wtmi1j1=dxmj/dx*dymi/dy;


            Tm(m)=T(i, j)*wtmij + T(i+1, j)*wtmi1j+T(i, j+1)*wtmij1+T(i+1, j+1)*wtmi1j1;
        end
    end
    if (nsteps>1)
        %         DT=T-T0;
        for m=1:1:nm
            i=fix(((ym(m)-yp(1)))/dy)+1;
            j=fix(((xm(m)-xp(1)))/dx)+1;

            dxmj=(xm(m)-xp(j));
            dymi=(ym(m)-yp(i));

            wtmij=(1-dxmj/dx)*(1-dymi/dy);
            wtmi1j=(1-dxmj/dx)*dymi/dy;
            wtmij1=dxmj/dx*(1-dymi/dy);
            wtmi1j1=dxmj/dx*dymi/dy;

            Tm(m)=Tm(m)+(DT(i, j))*wtmij + (DT(i+1, j))*wtmi1j+(DT(i, j+1))*wtmij1+(DT(i+1, j+1))*wtmi1j1;
        end
    end

    timesum=timesum+dt;

    figure(1); clf
    pcolor(xvy, yvy, RHO);
    title('density RHO')
    colormap('Jet')
    shading flat
    colorbar
    axis ij

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

    figure(7); clf
    pcolor(xp, yp, T) ; colormap('Jet')
    title('temperature T')
    shading flat
    colorbar
    axis ij

    figure(8); clf
    pcolor(xp, yp, T0) ; colormap('Jet')
    title('temperature T0')
    shading flat
    colorbar
    axis ij

    disp(p(27, 12));
    disp(vx(27, 12));
    disp(vy(27, 12));
    disp(T(27, 12));
    disp(dt);
    disp(timesum);

end

