% Solving of 1D Poisson eqn d2FI/dx2=2x^2-0.5*x+exp(x)
% with finite differences
% on a regular grid

% 0 Clearing variables and figures
clear all % clearing memory
clf % clearing figures

% 1 Define numerical model
xsize=1; % Horizontal model size, m
Nx=101; % Number of grid points
dx=xsize/(Nx-1);% Grid step, m,this is delta x
x=0:dx:xsize; % Coordinates of grid points, m

% 2 Define sparse matrix L() and vector R()
L=sparse(Nx, Nx); % number of equations, then number of unknowns, coefficients in the left hand side
R=zeros(Nx,1); % Lines, columns, values in the right hand side

% 3 Compose matrices L(), R()
% Going through all points of the grid
for j=1:1:Nx
    % Discriminate between BC-points and internal points
    if (j==1 || j==Nx)
        % BC equation 1*FI(j)=0, so for the first and last row
        L(j, j)=1; % Coefficient
        R(j, 1)=0; % Right hand side
    else
        % Poisson eqn d2FI/dx2=1
        %    FI(J-1)          FI(j)               FI(j+1)
        % ---[j-1]-------------[j]-----------------[j+1]
        %
        % d2FI/dx2(j)=(FI(j-1)-2*FI(j)+FI(j+1))/dx^2=2x^2-0.5*x+exp(x)
        % Left hand side
        % here I put the right hand side directly in the coefficient part
        % devided by 2x^2-x/2+exp(x)
        % x in the formula = (j-1)dx
        % !!! could be more simple
%         L(j, j-1)=1/((2.*(j-1)^2*dx^2-0.5*(j-1)*dx+exp((j-1)*dx))*dx^2); % FI(j-1)
%         L(j, j)=-2/((2.*(j-1)^2*dx^2-0.5*(j-1)*dx+exp((j-1)*dx))*dx^2); % FI(j)
%         L(j, j+1)=1/((2.*(j-1)^2*dx^2-0.5*(j-1)*dx+exp((j-1)*dx))*dx^2); % FI(j+1)
%         % Right hand side
%         R(j, 1)=1;
        L(j, j-1)=1/dx^2; % FI(j-1)
        L(j, j)=-2/dx^2; % FI(j)
        L(j, j+1)=1/dx^2; % FI(j+1)
        % Right hand side
        R(j, 1)=2*x(j)^2-0.5*x(j)+exp(x(j));
    end
end
% 3 Solve matrices
FI=L\R;
% 4 Visualise solution
figure(1)
plot(x,FI,'-o r')
% Compare with analytical solution
Nxa=1001; % number of points
dxa=xsize/(Nxa-1); % step between points
xa=0:dxa:xsize; % coordinates
FIa=exp(xa)+(xa.^4)/6-(xa.^3)/12+((1-exp(xsize))/xsize+(xsize^2-2.*xsize^3)/12)*xa-1;
hold on % continue plotting
plot(xa, FIa, '- b')
aaa(1,1)=FI(12);% -0.0819838745668572 numerical
aaa(2,1)=FIa(120);% -0.0881292935614402 analytical

