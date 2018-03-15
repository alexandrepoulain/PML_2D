clear; clc; close all;
%% PML 2D enhanced (2nd order) bar Test
% Parameters of input Ricker Wave
tp = 3; % Fundamental period
ts = 3; % Time shift
A = 1; % Amplitude

% Material parameters
rho_1 = 1; % density
mu_1 = 1; % 1st Lamé Coefficient 
nu_1 = 0.25; % Poisson modulus
E_1 = mu_1*(2*(1+nu_1)); % Young modulus

cs = (mu_1/rho_1)^0.5 ; % Shear wave velocity
lmbd_s = cs*tp; % wave length of shear wave

L1 = 250; % Size of the soil
h = 5;          % width of the plane
Lef = 5;    % Length of finite element
Nef = floor(L1 / Lef) ; % number of finite element in SDA
Nex = Nef; % nb elem in x direction
Ney = 1;  % nb elem in y direction
% Start coordinates
% down left
cx1 = 0; 
cy1 = 0;
% Up right
cx2 = L1; 
cy2 = h;
% Vector of coordinates in x and y direction
xp = cx1:(cx2-cx1)/Nex:cx2;
yp = cy1:(cy2-cy1)/Ney:cy2;
% Coordinate matrix
ind = 1;
for ii = 1:(Ney+1)
    for jj = 1:(Nex+1)
        XA(ind,:) = [xp(1,jj) yp(1,ii)];
        ind  = ind+1;
    end
end
% Connectivity
ip = 1 ; ll = 1 ;
for ii = 1:Ney
    for jj = 0:(Nex-1)
        TA(ll,:) = [(ip+jj) (ip+jj+1) (ip+jj+Nex+2) (ip+jj+Nex+1) 1] ;
        ll = ll + 1 ;
    end
    ip = ip + (Nex+1) ;
end

%% Space discretization (4-nodes bilinear elements)
Le_1 = 5; % Length of element 
dof = 2; % Per element
Nx_1 = Nef+1; Ny_1 = Nef+1; % number of nodes in x and y directions
% Square elements
a = Lef;
b = Lef;
% Nodal shape functions
N1 = @(xi, eta) 1/4*(1-xi)*(1-eta);
N2 = @(xi, eta) 1/4*(1+xi)*(1-eta);
N3 = @(xi, eta) 1/4*(1+xi)*(1+eta);
N4 = @(xi, eta) 1/4*(1-xi)*(1+eta);

C = E_1/(1-nu_1) * [1 nu_1 0;
                    nu_1 1 0;
                    0 0 (1-nu_1)/2];

N = @(xi,eta) [N1(xi,eta) 0 N2(xi,eta) 0 N3(xi,eta) 0 N4(xi, eta) 0;
    0 N1(xi,eta) 0 N2(xi,eta) 0 N3(xi, eta) 0 N4(xi,eta)];

% Strain matrix
B = @(xi,eta) [-(1-eta)/4 0 (1-eta)/4 0 (1+eta)/4 0 -(1+eta)/4 0 ;
    0 -(1-xi)/4 0 -(1+xi)/4 0 (1+xi)/4 0 (1-xi)/4 ;
    -(1-xi)/4 -(1-eta)/4 -(1+xi)/4 (1-eta)/4 (1+xi)/4 (1+eta)/4 (1-xi)/4 -(1+eta)/4];
% Derivative of N
dN = @(xi,eta)[-(1-xi)/4 -(1-eta)/4 -(1+xi)/4 (1-eta)/4;
    (1+xi)/4 (1+eta)/4 (1-xi)/4 -(1+eta)/4];

%% Gauss Quadrature 2 points in each direction
xi1 = -1/sqrt(3);
xi2 = 1/sqrt(3);
eta1 = -1/sqrt(3);
eta2 = 1/sqrt(3);
w1 = 1;
w2 = 1;
% Matrices construction
% mass 
M = assembleM(XA,TA,dof,Nef,rho_1);

% stiffness
K = zeros((Nef+1)*dof, (Nef+1)*dof);
for i = 1:Nef
   % Position of the nodes
   Xe = XA(TA(i,1:4),:);
   Te = TA(i,:);
   for k=1:4
       
       ke = (w1*(w1 *(B(xi1,eta1))' * C * B(xi1,eta1) + w2 *(B(xi1,eta2))'* C*B(xi1,eta2)) ...
           + w2*(w1 *(B(xi2,eta1))' * C * B(xi2,eta1) + w2 *(B(xi2,eta2))'* C*B(xi2,eta2))); 
   end
   K(i*dof:(i+1)*dof-1,i*dof:(i+1)*dof-1) = ke;
end

%% Time evolution
% parameters
dt_CFL = min(a,b) / sqrt(E_1/rho_1) ;
dt_1 = dt_CFL * 0.8;  
t_fin = 8.;
ntps_1 = round(t_fin/dt_1);
gamma1=0.5; beta1=0; % Explicit Newmark
%% Initializing
u1 = zeros(Nex+1,Ney+1);
u1_t = zeros(Nex+1,Ney+1);
u1_tt = zeros(Nex+1,Ney+1);

%% Plot
[X_plot,Y_plot] = meshgrid(0:L1/Nex:L1,0:h/Ney:h);
surf(X_plot,Y_plot,zeros(size(X_plot)));
shading interp;colormap('jet');
set(gcf, 'Position', get(0, 'Screensize'));
axis([0 L1 0 h -1 1]);
axis(gca,'on');
F(ntps_1) = struct('cdata',[],'colormap',[]);
F(1) = getframe(gcf);
for tt = 2:ntps_1
    u1 = 0;
end






