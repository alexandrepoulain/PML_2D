clear; clc; close all;
%% PML 2D enhanced (2nd order) bar Test
% Parameters of input Ricker Wave
tp = 3; % Fundamental period
ts = 3; % Time shift
A = 1; % Amplitude
b_scal = 1;

%% Material parameters
% SDA
rho_1 = 1700; % density
nu_1 = 0.24; % Poisson modulus
E_1 = 1e07; % Young modulus
mu_1 = E_1/(2*(1+nu_1)); % 1st Lamé Coefficient
lmbd_1 = E_1*nu_1/((1+nu_1)*(1-2*nu_1)); % 2nd Lamé Coefficient
k_1 = lmbd_1+2/3*mu_1;   % Bulk modulus
cs = (mu_1/rho_1)^0.5 ; % Shear wave velocity
cp = sqrt((lmbd_1+2*mu_1)/rho_1); % Pressure wave velocity
lmbd_s = cs*tp; % wave length of shear wave
L1 = 250; % Size of the soil
h = 5; % width of the plane
Lef = 5; % Length of finite element
Nef = floor(L1 / Lef) ; % number of finite element in SDA
Nex = Nef; % nb elem in x direction
Ney = 1; % nb elem in y direction
% Coordinates and connectivity matrix
[XA,TA] = topology(0,0,L1,h,Nex,Ney);
% Space discretization (4-nodes bilinear elements)
Le_1 = 5; % Length of element 
dof = 2; % Per element
% Square elements
a = Lef;
b = Lef;
% SDB
L2 = 200;   % Length of the PML
Le_2 = 5;   % Length of an element
NEF_2 = floor( L2 / Le_2 ); % Number of elements
% Number of nodes
NSx_2 = NEF_2 ; % number of nodes in x-direction
NSy_2 = 1 ; % number of nodes in y-direction
%% Time evolution
% parameters
dt_CFL = min(a,b) / sqrt(E_1/rho_1) ;
dt1 = dt_CFL*0.8;
dt2 = dt_CFL*0.8;
t_fin = 100.;
ntps_1 = round(t_fin/dt1);
gamma1=0.5; beta1=0; % Explicit Newmark
gamma2=0.5; beta2=0.25; % Implicit Newmark
%% Matrices construction
%%%%%%%%%%%%%%% SDA %%%%%%%%%%%%%%%
% mass 
MA = assembleM(XA,TA,dof,Nef,rho_1);
% stiffness
KA = assembleK(XA,TA,dof,E_1,nu_1,Nef);
% Lumped mass matrix 
% nombre de ddls total
NA = dof * size(XA,1) ;
diagM = zeros(NA,1) ;
for i=1:NA 
  diagM(i) = sum(MA(i,1:NA)) ;
end
%%%%%%%%%%%%%%% SDB %%%%%%%%%%%%%%%
n = 2; % order for attenuation functions
Rpp = 1e-3; % reflexion wanted
alpha_kucu = ((cp/cs)*(n+1))/(2*L2)*log(1/Rpp); % attenuation coefficient calculated by Kucukcoban
a0 = [alpha_kucu,alpha_kucu]; % coefficient of attenuation evanescent waves
b0 = [alpha_kucu,alpha_kucu]; % coefficient of attenuation propagating waves
x0 = L1; % start of the PML
% Coordinates and connectivity matrix
[XB,TB] = topology(L1,0,L1+L2,h,NSx_2,NSy_2);
% Matrices of PML
% mass
MB = assembleM_PML(rho_1,XB,TB,a0,x0,L2,h,n,dof);
% stiffness
KB = assembleK_PML(rho_1,cs,b_scal,XB,TB,a0,b0,x0,L2,h,n,dof);
% damping
CB = assembleC_PML(rho_1,cs,b_scal,XB,TB,a0,b0,x0,L2,h,n,dof);
% effective stiffness
KKB = assemble_effKB_PML(k_1,mu_1,cs,b_scal,dt2,XB,TB,a0,b0,x0,L2,h,n,dof);
% effective damping
CCB = assemble_effCB_PML(k_1,mu_1,cs,b_scal,dt2,XB,TB,a0,b0,x0,L2,h,n,dof);
%% Nodes at the interface
[interfaceA, interfaceB] = count_nodes_interface(L1,XA,XB);


%% Initializing
nnA = size(XA,1);
u1 = zeros(nnA*dof,1);
v1 = zeros(nnA*dof,1);
a1 = zeros(nnA*dof,1);
%% Nodes with load
load_nodes = zeros(1,2);
compt = 1;
for i=1:nnA
  if (XA(i,1) == 0)
      load_nodes(compt) = i;
      compt = compt+1;
  end
end
Forn1  = zeros(dof*nnA,1) ; 
for i=1:length(load_nodes)
  Forn1((load_nodes(i)-1)*dof+1)=1 ;
end
%% Boundary conditions
% On SDA: fixed end
BCAx = [51 102];
BCAy = [51 102];
%% Plot
% Reshape for plotting
subplot(2,2,1);
coord1=0:Le_1:L1;
axis([0 L1 -5 5]);
set(gca,'NextPlot','replacechildren');
size1= size(coord1,2);
plot(coord1,u1(1:2:102),'k-','linewidth',2);

subplot(2,2,2);
coord1=0:Le_1:L1;
axis([0 L1 -5 5]);
set(gca,'NextPlot','replacechildren');
plot(coord1,u1(2:2:102),'k-','linewidth',2);

subplot(2,2,[3,4]);
X_plot(1:51,1) = XA(1:51,1); X_plot(1:51,2) = XA(52:102,1);
Y_plot(1:51,1) = XA(1:51,2); Y_plot(1:51,2) = XA(52:102,2);
surf(X_plot,Y_plot,zeros(size(X_plot)));
axis([0 L1 0 h -1 1]);
axis(gca,'on');
F(ntps_1) = struct('cdata',[],'colormap',[]);
F(1) = getframe(gcf);

t1=0;
for tt = 2:ntps_1
    % Ricker force
    f1 = F_ricker(1E6,t1); % force
    Fext1 = f1*Forn1;      % application
    % Predictors
    up1 = u1+(dt1*v1)+(((0.5-beta1)*(dt1^2))*a1);  
    vp1 = v1+(((1.0-gamma1)*dt1)*a1);
    % Righthand term
    rhs = Fext1-(KA*up1); 
    % acceleration
    a1 = rhs./diagM;
    % update
    v1  = vp1+((gamma1*dt1)*a1); 
    u1  = up1+((beta1*(dt1^2))*a1); 
    
    
    
    %% Plot
    subplot(2,2,1);
    plot(coord1,u1(1:2:102),'k-','linewidth',2);
    subplot(2,2,2);
    plot(coord1,u1(2:2:102),'k-','linewidth',2);
    subplot(2,2,[3,4]);
    XA_new(:,1) = XA(:,1) + u1(1:2:204);
    XA_new(:,2) = XA(:,2) + u1(2:2:204);
    Z(1:51,1) = u1(1:2:102); Z(1:51,2) = u1(103:2:204);
    % Reshape for plotting
    X_plot(1:51,1) = XA_new(1:51,1); X_plot(1:51,2) = XA_new(52:102,1);
    Y_plot(1:51,1) = XA_new(1:51,2); Y_plot(1:51,2) = XA_new(52:102,2);
    surf(X_plot,Y_plot,zeros(size(X_plot)),Z);
    F(tt)=getframe;
    % time step
    t1=t1+dt1;
    tt
end







