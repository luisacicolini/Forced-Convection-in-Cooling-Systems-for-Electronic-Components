%%%% DATA

boundary.south = 'Dirichlet';
boundary.north = 'Robin';
boundary.east = 'Dirichlet';
boundary.west = 'Dirichlet';

% Values for Dirichlet BC

TD.north = 298;
TD.south= 298;
TD.west= 298;
TD.east= 298;


%% Gemetrical Parameters %%

l = 2;                 % Length of the domain in x 
% h = 2;                 % Length of the domain in y
f = 2;                 % Length of the domain in z
dimX = 20;            % Number of nodes in x
% dimY = 3;            % Number of nodes in y
z1 = 20;                % Number of nodes in z (board)
z2 = 20;                % Number of nodes in z (fluid)
dimZ = z1+z2;          % Number of nodes in z (total)
% dimXY= dimX*dimY;
dimXZ= dimX*dimZ;
dimXz1= dimX*z1;
dimXz2= dimX*z2;

%% Data %%
rho= [1.850 1.225]';               % Density of the [board fluid]
c= [0.396 1]';               % Heat capacity [board fluid]
D= [1.446 0.08]';               % Diffusivity [board fluid]
%% Velocity of convection
m= 2;
v_x= zeros(z2,1);             % Velocity of forced convection
v_x= v_x + m*(1:z2)';             % Velocity of forced convection
% % v_x(8:end)= 0;

% Values of Neumann BC
beta = 1;                % Heat flux
% Values of Robin BC
alpha = 20;               % Convective heat transfer coefficient
Tinf = 298;               % Temperature of the surrounding fluid

% Thermal conductivity Coefficient
heat_conduc='homogeneous'  % Choose: 1) homogeneous, 2) non homogeneous
% Kval = [rho(1)*c(1)*D(1) rho(2)*c(2)*D(2)];  % Thermal conductivity (homogeneous region)
% Kval =[0.2e3 24.e-3];
Kval =[1.059 0.02509];
% KnH = 10;              % Thermal conductivity for the non-homogeneus region
% yk = [0.5 0.8];   % vertical coordinates of non-homogeneous region
% xk = [1 1.5];     % horizontal coordinates of non-homogeneous region

%% Source parameters. Circle ==>  (x-c)^2+(z-d)^2<= R
source= 'yes' %Choose:   1) yes   2) no

% Board

w1 = 333;      w2 = 333;
c1 = 0.65;       c2 = 1.35;         %[0 l]
d1 = 0.9*f/2;   d2 = d1;    %[0 f*] inside the board
R1 = 0.15;       R2 = 0.15;




% % % Fluid
% % 
% % w3 = 398;
% % c3 = l/2;
% % d3 = 0.75*f;
% % R3 = 0.15;

% Temperature first values
Tmean= (TD.west+TD.east+TD.north+TD.south)/4;
Tnb = Tmean*ones(z1, dimX);
Tnf = Tmean*ones(z2, dimX);
Tnb(:,1) = TD.west*strcmp(boundary.west,'Dirichlet')+...
    (1-strcmp(boundary.west,'Dirichlet'))*Tnb(:,1);
Tnb(:,end) = TD.east*strcmp(boundary.east,'Dirichlet')+...
    (1-strcmp(boundary.east,'Dirichlet'))*Tnb(:,end);
Tnb(1,:) = TD.north*strcmp(boundary.north,'Dirichlet')+...
    (1-strcmp(boundary.north,'Dirichlet'))*Tnb(1,:);
Tnb(end,:) = TD.south*strcmp(boundary.south,'Dirichlet')+...
    (1-strcmp(boundary.south,'Dirichlet'))*Tnb(end,:);
