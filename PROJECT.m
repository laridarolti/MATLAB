% 0 clear memory and figures
clear all % memory
clf % figures


% DEFINE NUMERICAL MODEL
xsize = 500000; ysize = 400000; % horizontal, vertical, m
nx=101;ny=81; % grid resolution in horizontal, vertical direction
dx=xsize/(nx-1); dy=ysize/(ny-1); % horizontal, vertical grid step in m
x=0:dx:xsize+dx; y=0:dy:ysize+dy; % coordinates of grid points
xp=-dx/2:dx:xsize+dx/2; yp=-dy/2:dy:ysize+dy/2;
xvx=0:dx:xsize+dx; yvx=-dy/2:dy:ysize+dy/2;
xvy=-dx/2:dx:xsize+dx/2; yvy=0:dy:ysize+dy;
gy=10; % defining gravitational constant
r=20000; % defining radius of plume
BC_type=-1; % free slip=-1, no slip =1
RHOmelt=2500;
ETA_melt=10;
RHO = zeros(ny+1, nx+1);RHOCP = zeros(ny+1, nx+1);ETA_P = zeros(ny+1, nx+1);ETA_B = zeros(ny+1, nx+1);kx= zeros(ny+1, nx+1);ky= zeros(ny+1, nx+1);T0= zeros(ny+1, nx+1);alpha=zeros(ny+1, nx+1);HR=zeros(ny+1, nx+1);dt_rate=zeros(ny+1, nx+1);DT_subgrid=zeros(ny+1, nx+1);
NT=(nx+1)*(ny+1)*1; % total number of unknowns
LT=sparse(NT, NT); % coeff in the left side of equations
RT=zeros(NT,1); % values for the right hand side of eqns
dt=6e+10; dispmax=0.5; DTmax=50;%dtthermal= min(dx,dy)^2/(4*max([4/3350/800 3/3300/1000  2/3200/1100])); % CHECK dtthermal
timesum=0; ETAmin=1e+18; ETAmax=1e+24;

% define marker arrays
nxm=(nx-1)*4; nym=(ny-1)*4;nm=nxm*nym;% number of markers in horizontal direction
xm=zeros(nm, 1); ym=zeros(nm, 1);% x coordinates of markers
dxm=xsize/nxm; dym=ysize/nym;% average horizontal distance between markers in m
RHO0m=zeros(nm, 1);RHOm=zeros(nm, 1);RHOinitialm=zeros(nm, 1);% density value for markers in kg/m3
ETA0m=zeros(nm, 1);ETAm=zeros(nm, 1);ETAinitialm=zeros(nm, 1); % viscosity values for markers
CPm=zeros(nm, 1);RHOCPm=zeros(nm, 1);km=zeros(nm, 1);Tm=zeros(nm, 1); a=13/600; b=4173;
T0m=zeros(nm, 1);dt_ratem=zeros(nm, 1);dT_relaxedm=zeros(nm, 1);dT_subgridm=zeros(nm, 1);dTm=zeros(nm, 1);
HRm=zeros(nm, 1);alpham=zeros(nm, 1);bettam=zeros(nm, 1);strengthm=zeros(nm, 1);pm=zeros(1,nm);
EIIm=zeros(1,nm);
% for wet/ dry olivine
Eam=zeros(nm, 1);Vam=zeros(nm, 1);power_nm=zeros(nm, 1);Adm=zeros(nm, 1);

m=1;
for jm=1:1:nxm
    for im=1:1:nym
        % coordinates
        xm(m)=dxm/2+(jm-1)*dxm +(rand-0.5)*dxm; %horizontal
        ym(m)=dym/2+(im-1)*dym +(rand-0.5)*dym; %vertical

        % I
        if (ym(m)<50000)
            RHOinitialm(m) =1;HRm(m)=0;km(m)=3000;ETA0m(m) = 1e+18;Tm(m)=273;CPm(m)=3300000;%RHO0m(m) =1;

            % II
        elseif (ym(m)>=50000 && ym(m)<100000 )
            RHOinitialm(m) =3400;RHO0m(m) =3400;alpham(m)=3e-5;bettam(m)=1e-11;HRm(m)=2e-8;
            Tm(m)=-1027+ym(m)*13/500; % increases gradually
            if (ym(m)>(xm(m)-180000) && ym(m)<(xm(m)-120000 ) && ym(m)>(-5*xm(m)+1200000))
                % y=-5x+12*10^5 is the third line in that triangle i.e. the line connecting the two parallel lines
                c=ym(m)-xm(m);Tm(m)=a*c+b;
            end
            if (Tm(m)>-77)
                km(m)=0.73+1293/(Tm(m)+77);
            end
            CPm(m)=1000;strengthm(m)=1e+8;ETA0m(m) = 1e+23;
            Eam(m)=532000;Vam(m)=8e-6;power_nm(m)=3.5;Adm(m)=2e-17;% DRY OLIVINE


            % III
        elseif (ym(m)>=100000 && ym(m)< 150000 && ym(m)>(xm(m)-180000) && ym(m)<(xm(m)-120000 ))
            RHOinitialm(m) =3400;RHO0m(m) =3400;alpham(m)=3e-5;bettam(m)=1e-11;HRm(m)=2e-8;ETA0m(m) = 1e+23;CPm(m)=1000;strengthm(m)=2e+7;
            c=(ym(m)-xm(m));Tm(m)=a*c+b;
            if (Tm(m)>-77)
                km(m)=0.73+1293/(Tm(m)+77);
            end
            Eam(m)=532000;Vam(m)=8e-6;power_nm(m)=3.5;Adm(m)=2e-17;% DRY OLIVINE

            % IV
        elseif (ym(m)>=150000 && ym(m)< 250000 && ym(m)>(xm(m)-180000) && ym(m)<(xm(m)-120000 ))
            RHOinitialm(m) =3400;RHO0m(m) =3400;alpham(m)=3e-5;bettam(m)=1e-11;HRm(m)=2e-8;ETA0m(m) = 1e+23;CPm(m)=1000;strengthm(m)=1e+8;
            c=ym(m)-xm(m);Tm(m)=a*c+b;
            if (Tm(m)>-77)
                km(m)=0.73+1293/(Tm(m)+77);
            end
            Eam(m)=532000;Vam(m)=8e-6;power_nm(m)=3.5;Adm(m)=2e-17; % DRY OLIVINE

        else
            RHO0m(m) =3300;RHOinitialm(m) =3250;alpham(m)=3e-5;bettam(m)=1e-11;HRm(m)=3e-8;ETA0m(m) = 1e+20;CPm(m)=1000;strengthm(m)=5e+7;
            Tm(m)=1573;
            if (Tm(m)>-77)
                km(m)=0.73+1293/(Tm(m)+77);
            end
            Eam(m)=471000;Vam(m)=4e-6;power_nm(m)=4;Adm(m)=2e-21; % WET OLIVINE
        end
        RHOm(m)=RHOinitialm(m);RHOCPm(m)=RHOm(m)*CPm(m);ETAm(m)=ETA0m(m);
        m=m+1;
    end
end
aaa=zeros(nm, 1); % finding dt thermal = min(dx, dy)^2/4*K
for m=1:1:nm
    if (RHOCPm(m)~=0)
        aaa(m)=km(m)/RHOCPm(m);
    end
end
dtthermal= min(dx,dy)^2/(4*max(aaa));
nsteps=10;

for timestep=1:1:nsteps
    disp('timestep is'); disp(timestep);
    if (timestep>1)

        % INTERPOLATING P(to calculate rho), EII
        for m=1:1:nm
            i=fix(((ym(m)-yp(1)))/dy)+1;
            j=fix(((xm(m)-xp(1)))/dx)+1;
            dxmj=(xm(m)-xp(j));
            dymi=(ym(m)-yp(i));
            wtmij=(1-dxmj/dx)*(1-dymi/dy);
            wtmi1j=(1-dxmj/dx)*dymi/dy;
            wtmij1=dxmj/dx*(1-dymi/dy);
            wtmi1j1=dxmj/dx*dymi/dy;

            pm(m)=p(i, j)*wtmij + p(i+1, j)*wtmi1j+p(i, j+1)*wtmij1+p(i+1, j+1)*wtmi1j1;
            EIIm(m)=EII(i, j)*wtmij + EII(i+1, j)*wtmi1j+EII(i, j+1)*wtmij1+EII(i+1, j+1)*wtmi1j1;

            % CALCULATING RHOm AND ETAm
            if (ym(m)>=50000) % quantities stay unchanged for sticky air
                if(1+alpham(m)*(Tm(m)-273)~=0)
                    RHOm(m)=RHO0m(m)*(1+bettam(m)*(pm(m)-1e+5))/(1+alpham(m)*(Tm(m)-273));
                end
                RHOCPm(m)=RHOm(m)*CPm(m);
                if (Tm(m)~=0 && Adm(m)~=0 && power_nm(m)~=0)
                    ETAm(m)=0.5/Adm(m)^(1/power_nm(m))*EIIm(m)^(1/power_nm(m)-1)*exp((Eam(m)+pm(m)*Vam(m))/8.314/Tm(m)/power_nm(m));
                end
                if (ETAm(m)>ETAmax)
                    ETAm(m)=ETAmax;
                end
                if (EIIm(m)~=0 && ETAm(m)>0.5*strengthm(m)/EIIm(m))
                    ETAm(m)=0.5*strengthm(m)/EIIm(m); % plastic viscosity limit
                end
                if (ETAm(m)<ETAmin)
                    ETAm(m)=ETAmin;
                end
                if (Tm(m)>-77)
                    km(m)=0.73+1293/(Tm(m)+77);
                end
            end
        end
    end

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
    HR_WT=zeros(ny+1, nx+1);
    HR_SUM=zeros(ny+1, nx+1);
    alpha_WT=zeros(ny+1, nx+1);
    alpha_SUM=zeros(ny+1, nx+1);

    % calculate ETA_P, ETA_B, RHO (or RHOvy in vy points)
    for m=1:1:nm
        % I. ETA_P, Tdt and RHOCP, HR
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
        T_SUM(i+1, j)=T_SUM(i+1, j)+wtmi1j*Tm(m)*RHOCPm(m);
        T_SUM(i, j+1)=T_SUM(i, j+1)+wtmij1*Tm(m)*RHOCPm(m);
        T_SUM(i+1, j+1)=T_SUM(i+1, j+1)+wtmi1j1*Tm(m)*RHOCPm(m);
        %  calculating numerator and denominator for alpha
        alpha_SUM(i, j)= alpha_SUM(i, j)+wtmij* alpham(m);
        alpha_WT(i, j)= alpha_WT(i, j)+wtmij;
        alpha_SUM(i+1, j)= alpha_SUM(i+1, j)+wtmi1j* alpham(m);
        alpha_WT(i+1, j)= alpha_WT(i+1, j)+wtmi1j;
        alpha_SUM(i, j+1)= alpha_SUM(i, j+1)+wtmij1* alpham(m);
        alpha_WT(i, j+1)= alpha_WT(i, j+1)+wtmij1;
        alpha_SUM(i+1, j+1)= alpha_SUM(i+1, j+1)+wtmi1j1* alpham(m);
        alpha_WT(i+1, j+1)= alpha_WT(i+1, j+1)+wtmi1j1;
        % calculating numerator and denominator for HR
        HR_SUM(i, j)=HR_SUM(i, j)+wtmij*HRm(m);
        HR_WT(i, j)=HR_WT(i, j)+wtmij;
        HR_SUM(i+1, j)=HR_SUM(i+1, j)+wtmi1j*HRm(m);
        HR_WT(i+1, j)=HR_WT(i+1, j)+wtmi1j;
        HR_SUM(i, j+1)=HR_SUM(i, j+1)+wtmij1*HRm(m);
        HR_WT(i, j+1)=HR_WT(i, j+1)+wtmij1;
        HR_SUM(i+1, j+1)=HR_SUM(i+1, j+1)+wtmi1j1*HRm(m);
        HR_WT(i+1, j+1)=HR_WT(i+1, j+1)+wtmi1j1;

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
        % III. RHOvy, ky, alpha
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
        % weights for vx pts, use same variables, previous values not needed
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
            if (RHOCP_SUM(i, j)>0)
                T0(i, j)=T_SUM(i, j)/RHOCP_SUM(i, j);
            else
                T0(i, j)=T0(i, j);
            end
            if (HR_WT(i, j)>0)
                HR(i, j)=HR_SUM(i, j)/HR_WT(i, j);
            else
                HR(i, j)=HR(i, j);
            end
            if (alpha_WT(i, j)>0)
                alpha(i, j)=alpha_SUM(i, j)/alpha_WT(i, j);
            else
                alpha(i, j)=alpha(i, j);
            end
        end
    end

    % 2 defining global matrices L and R

    N=(nx+1)*(ny+1)*3; % total number of unknowns
    L=sparse(N, N); % coeff in the left side of equations
    R=zeros(N,1); % values for the right hand side of eqns
    dt=dt*1.2; % increase by 20% if timestep is too slow
    niter=3; % 2 iter per timestep
    Tdt=T0;
    for iter=1:1:niter
        for j=1:1:nx+1
            for i=1:1:ny+1
                gp=((j-1)*(ny+1)+i-1)*3+1;gvx=gp+1;gvy=gp+2;

                % continuity equation
                if (i==1 || j==1 || j==nx+1 || i==ny+1)
                    L(gp, gp)=1;
                    R(gp, 1)=0;
                elseif (i==2 && j==2)
                    L(gp, gp)=1;
                    R(gp, 1)=0; %1e+5; %3.3E9;
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

        S=L\R; % 4 solve matrices

        % 5 reload solution into vx, vy, p % define empty geometrical arrays
        p=zeros(ny+1, nx+1);vx=zeros(ny+1, nx+1);vy=zeros(ny+1, nx+1);
        for j=1:1:nx+1 % from 1 with step 1
            for i=1:1:ny+1
                gp=((j-1)*(ny+1)+i-1)*3+1;gvx=gp+1;gvy=gp+2;% define global index of the current equation S
                p(i, j)= S(gp); vx(i, j)= S(gvx);vy(i, j)= S(gvy);% reload S(gp) => p, vx, vy
            end
        end

        EXY=zeros(ny+1, nx+1); SXY=zeros(ny+1, nx+1);
        for j=1:1:nx % from 1 with step 1
            for i=1:1:ny
                EXY(i, j)=((vx(i+1, j)-vx(i, j))/dy+(vy(i, j+1)-vy(i, j))/dx)/2;SXY(i, j)=2*ETA_B(i, j)*EXY(i, j);
            end
        end

        EXX=zeros(ny+1, nx+1); SXX=zeros(ny+1, nx+1);EII=zeros(ny+1, nx+1);
        for j=2:1:nx % from 1 with step 1
            for i=2:1:ny
                EXX(i, j)=(vx(i, j)-vx(i, j-1))/dx;SXX(i, j)=2*ETA_P(i, j)*EXX(i, j);EII(i, j)=sqrt(EXX(i, j)^2+EXY(i, j)^2);
            end
        end

        HS=zeros(ny+1, nx+1);
        for j=2:1:nx
            for i=2:1:ny
                HS(i, j)=2*SXX(i, j)*EXX(i, j)+2*(EXY(i, j)*SXY(i, j)+EXY(i-1, j)*SXY(i-1, j)+EXY(i, j-1)*SXY(i, j-1)+EXY(i-1, j-1)*SXY(i-1, j-1))/4;
            end
        end

        HA=zeros(ny+1, nx+1);
        T=(T0+Tdt)/2;
        for j=1:1:nx+1
            for i=2:1:ny-1
                HA(i, j)=alpha(i, j)*T(i, j)*(RHO(i, j)*vy(i, j)+RHO(i-1, j)*vy(i-1, j))*gy/2;
            end
        end

        vxmax=max(max(abs(vx)));vymax=max(max(abs(vy)));
        if (dt*vxmax>dx*dispmax && iter<niter) %
            dt=dx*dispmax/vxmax;
        end
        if (dt*vymax>dy*dispmax && iter<niter ) %
            dt=dy*dispmax/vymax;
        end
        dt=min(dt, dtthermal);


        % 3 composing TEMP matrices
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
                    RT(gp, 1)=RHOCP(i, j)*T0(i, j)+dt*HR(i, j)+dt*HS(i, j)+dt*HA(i, j);
                end
            end
        end


        % 4 SOLVE MATRICES
        ST=LT\RT;
        % 5 RELOAD ST(gp) INTO Tdt(i, j)
        Tdt=zeros(ny+1, nx+1);
        for j=1:1:nx+1 % from 1 with step 1
            for i=1:1:ny+1
                gp=((j-1)*(ny+1)+i-1)*1+1;
                Tdt(i, j)= ST(gp);
            end
        end
        DT=Tdt-T0;
        if ( max(max(abs(DT)))>DTmax  && iter<niter)
            dt=dt*0.7*DTmax/max(max(abs(DT)));
        end
    end
    dt=min(dt, dtthermal);

    for j=2:1:nx+1
        for i=2:1:ny+1
            dt_rate(i, j)=exp(-dt*((ky(i, j)+ky(i-1, j))/dy^2+(kx(i, j-1)+kx(i, j))/dx^2)/RHOCP(i, j)); % CHECK THIS!!!
        end
    end

    % for interp dT subgrid from markers to nodal pts
    RHOCP_SUM=zeros(ny+1, nx+1);RHOCP_WT=zeros(ny+1, nx+1);T_SUM= zeros(ny+1, nx+1);

    % MOVING Tdt MARKERS
    if (timestep==1)
        for m=1:1:nm
            i=fix(((ym(m)-yp(1)))/dy)+1;
            j=fix(((xm(m)-xp(1)))/dx)+1;
            dxmj=(xm(m)-xp(j));
            dymi=(ym(m)-yp(i));
            wtmij=(1-dxmj/dx)*(1-dymi/dy);
            wtmi1j=(1-dxmj/dx)*dymi/dy;
            wtmij1=dxmj/dx*(1-dymi/dy);
            wtmi1j1=dxmj/dx*dymi/dy;
            Tm(m)=Tdt(i, j)*wtmij + Tdt(i+1, j)*wtmi1j+Tdt(i, j+1)*wtmij1+Tdt(i+1, j+1)*wtmi1j1;
        end

    else
        for m=1:1:nm
            i=fix(((ym(m)-yp(1)))/dy)+1;
            j=fix(((xm(m)-xp(1)))/dx)+1;
            dxmj=(xm(m)-xp(j));
            dymi=(ym(m)-yp(i));
            wtmij=(1-dxmj/dx)*(1-dymi/dy);
            wtmi1j=(1-dxmj/dx)*dymi/dy;
            wtmij1=dxmj/dx*(1-dymi/dy);
            wtmi1j1=dxmj/dx*dymi/dy;
            T0m(m)=(T0(i, j))*wtmij + (T0(i+1, j))*wtmi1j+(T0(i, j+1))*wtmij1+(T0(i+1, j+1))*wtmi1j1;
            dTm(m)=Tm(m)-T0m(m); % dT marker nodal = Tm - T0 marker nodal
            dt_ratem(m) = (dt_rate(i, j))*wtmij + (dt_rate(i+1, j))*wtmi1j+(dt_rate(i, j+1))*wtmij1+(dt_rate(i+1, j+1))*wtmi1j1; % dt_rate on markers
            dT_relaxedm(m)= dTm(m)*dt_ratem(m);

            % calculating dT subgrid + interp to nodal pts
            dT_subgridm(m)=dT_relaxedm(m)-dTm(m);

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
            T_SUM(i, j)=T_SUM(i, j)+wtmij*dT_subgridm(m)*RHOCPm(m);
            T_SUM(i+1, j)=T_SUM(i+1, j)+wtmi1j*dT_subgridm(m)*RHOCPm(m);
            T_SUM(i, j+1)=T_SUM(i, j+1)+wtmij1*dT_subgridm(m)*RHOCPm(m);
            T_SUM(i+1, j+1)=T_SUM(i+1, j+1)+wtmi1j1*dT_subgridm(m)*RHOCPm(m);

            %             Tm(m)=Tm(m)+dT_subgridm(m);

        end

        for j=1:1:nx+1
            for i=1:1:ny+1
                if (RHOCP_SUM(i, j)>0)
                    DT_subgrid(i, j)=T_SUM(i, j)/RHOCP_SUM(i, j);
                else
                    DT_subgrid(i, j)=DT_subgrid(i, j);
                end
                DT(i, j)=DT(i, j)-DT_subgrid(i, j);
            end
        end

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

    % INTERPOLATING VX TO VXM
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
            %             disp('j is');
            %             disp(j);
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
        xm(m)=xm(m)+vxm_eff(m)*dt;
        ym(m)=ym(m)+vym_eff(m)*dt;
    end


    timesum=timesum+dt;

    disp('p is'); disp(p(27, 12));
    disp('vx is'); disp(vx(27, 12));
    disp('vy is'); disp(vy(27, 12));
    disp('Tdt is'); disp(Tdt(27, 12));
    disp('dt is'); disp(dt);
    disp('timesum is'); disp(timesum);


    figure(1); clf
    pcolor(xvy, yvy, RHO);
    title('density RHO')
    colormap('Jet')
    shading flat
    colorbar
    axis ij

    figure(2); clf
    pcolor(xp, yp, log(ETA_P)) ; colormap('Jet')
    title('ETA_P viscosity in pressure points')
    shading flat
    colorbar
    axis ij

    figure(3); clf
    pcolor(x, y, log(ETA_B)) ; colormap('Jet')
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
    pcolor(xp, yp,T0) ; colormap('Jet')
    title('temperature T0')
    shading flat
    colorbar
    axis ij

    figure(8); clf
    pcolor(xp, yp, RHOCP) ; colormap('Jet')
    title('RHOCP')
    shading flat
    colorbar
    axis ij

    figure(9); clf
    pcolor(xp, yp, HA) ; colormap('Jet')
    title('HA')
    shading flat
    colorbar
    axis ij

    figure(10); clf
    pcolor(xp, yp, HR) ; colormap('Jet')
    title('HR')
    shading flat
    colorbar
    axis ij

    figure(11); clf
    pcolor(xp, yp, EII) ; colormap('Jet')
    title('EII')
    shading flat
    colorbar
    axis ij

    figure(12); clf
    pcolor(xp, yp, Tdt) ; colormap('Jet')
    title('Tdt')
    shading flat
    colorbar
    axis ij

    figure(13); clf
    pcolor(xp, yp, HS) ; colormap('Jet')
    title('HS')
    shading flat
    colorbar
    axis ij

    figure(14); clf
    pcolor(xvx, yvx, log(kx)) ; colormap('Jet')
    title('log(kx)')
    shading flat
    colorbar
    axis ij
end

% reference values
disp(p(7, 5));
disp(vx(7, 5));
disp(vy(7, 5));