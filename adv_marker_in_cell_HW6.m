% solving 1D advection eqn using marker in cell method
% dB/dt=-vx*dB/dx
% on a regular grid

% clear mem and fig
clear all
clf

% define Eulerian grid
xsize=100; % horizontal model size, m
nx=101; % model resolution
dx=xsize/(nx-1); % grid step size, m
x=0:dx:xsize; % coord of eulerian grid pts
dt=0.1; % timestep, s

% define Lagrangian markers
nxm=500; % model resolution
dxm=xsize/nxm; % grid step size, m
xm=dxm/2:dxm:xsize-dxm/2; % coord of eulerian grid pts

% define initial parameter Bm distribution for markers
Bm=zeros(1, nxm);
for m=1:1:nxm
    if (xm(m)>xsize*0.4 && xm(m)< xsize*0.6)
        Bm(m)=3300; % density wave
    else
        Bm(m)=3200; % bacground density
    end
end

% initialising B0, vx velocity and vxm velocity in markers
 B0=zeros(1, nx);
 vx=zeros(1, nx);
 vxm=zeros(1, nxm);

% time stepping
for timestep=1:1:100


% interpolation from markers - Bm to nodal points - B0  
for j=1:1:nx
  SUM_weight=0; % for every point j, the sum of weights initialised with 0
  SUM_Bm=0;     % for every j sum of weight times Bm initialised with 0
  for m=1:1:nxm  
        if (abs(xm(m)-x(j))<=dx)                   % check if is within dx of a point j, if yes, then add Bm contributes to B0(j)
            wtmj=(1-abs(xm(m)-x(j)))/dx;           % weight formula
            SUM_Bm=SUM_Bm+wtmj*Bm(m);
            SUM_weight=SUM_weight+wtmj;
        end
  end
    B0(j)=SUM_Bm/SUM_weight;
    vx(j)= 1 + sin(x(j)/xsize*pi*5)*0.1*(B0(j)-3250)/50;   % vx velocity in nodal points
end


% interpolationg vx to vxm
for m=1:1:nxm
    if (xm(m)>xsize)
        xm(m)=xm(m)-xsize;
    end
    if (xm(m)<0)
        xm(m)=xm(m)+xsize;
    end
    j=fix(xm(m)/dx)+1;
    vxm(m)=vx(j)*(1-abs(xm(m)-x(j)))/dx + vx(j+1)*abs(xm(m)-x(j))/dx; %delta xm(j)= xm-x for the point

    % moving markers
    xm(m)=xm(m)+vxm(m)*dt;
end

 
% figure(1)
% plot(x, B0,'-o r')
% axis([0 xsize 3150 3350])
% figure(2)
% plot(x, vx,'-o r')
% axis([0 xsize 0.9 1.1])

end

figure(1)
plot(x, B0,'-o r')
axis([0 xsize 3150 3350])
figure(2)
plot(x, vx,'-o r')
axis([0 xsize 0.9 1.1])


disp(B0(51));

disp(vx(51));