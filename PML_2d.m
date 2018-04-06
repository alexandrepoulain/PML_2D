clear; clc; close all;
%% PML 2D enhanced (2nd order) bar Test
% Parameters of input Ricker Wave
tp = 3; % Fundamental period
ts = 3; % Time shift
A = 1; % Amplitude


%% Material parameters
% SDA
rho_1 = 1700; % density
nu_1 = 0.24; % Poisson modulus
E_1 = 1e07; % Young modulus
mu_1 = E_1/(2*(1+nu_1)); % 1st Lam� Coefficient
lmbd_1 = E_1*nu_1/((1+nu_1)*(1-2*nu_1)); % 2nd Lam� Coefficient
k_1 = E_1/3/(1-2*nu_1);   % Bulk modulus
cs = (mu_1/rho_1)^0.5 ; % Shear wave velocity
cp = sqrt((lmbd_1+2*mu_1)/rho_1); % Pressure wave velocity
lmbd_s = cs*tp; % wave length of shear wave
b_scal = cs;
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
Nef_2 = floor( L2 / Le_2 ); % Number of elements
% Number of nodes
NSx_2 = Nef_2 ; % number of nodes in x-direction
NSy_2 = 1 ; % number of nodes in y-direction
%% Time evolution
% parameters
dt_CFL = min(a,b) / sqrt(E_1/rho_1) ;
dt1 = dt_CFL*0.8;
dt2 = dt_CFL*0.8;
t_fin = 100.;
ntps_1 = round(t_fin/dt1);
gamma1=0.5; beta1=0.25;    % Implicit Newmark
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
alpha_kucu = sqrt(E_1/rho_1)*(n+1)/(2*L2)*log(1/Rpp); % attenuation coefficient calculated by Kucukcoban
alpha_kucu = 0;
a0 = [alpha_kucu,alpha_kucu]; % coefficient of attenuation evanescent waves
b0 = [alpha_kucu,alpha_kucu]; % coefficient of attenuation propagating waves
x0 = L1; % start of the PML
% Coordinates and connectivity matrix
[XB,TB] = topology(L1,0,L1+L2,h,NSx_2,NSy_2);
% From local to global and in 2D
global_pos_elemB = zeros(4*dof,size(TB,1));
for i=1:size(TB,1)
    nnodes = size(TB,2)-1;  
    for s = 1:nnodes
        for j = 1:dof
          global_pos_elemB((s-1)*dof+j,i) = (TB(i,s)-1)*dof + j ;
        end
    end
end
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
% dynamical operator for SDB
M_eff= (MB+gamma2*dt2*(CB+CCB)+ beta2*(dt2)^2*(KB+KKB));
M_eff=sparse(M_eff);
% inverse
inv_M_eff=inv(M_eff) ;
%% Nodes at the interface
[interfaceA, interfaceB] = count_nodes_interface(L1,XA,XB);


%% Initializing
nnA = size(XA,1); % number of nodes SDA
u1 = zeros(nnA*dof,1); % displacement SDA
v1 = zeros(nnA*dof,1); % velocity SDA
a1 = zeros(nnA*dof,1); % acceleration SDA
nnB = size(XB,1); % number of nodes SDB
u2 = zeros(nnB*dof,1); % displacement SDB
v2 = zeros(nnB*dof,1); % velocity SDB
a2 = zeros(nnB*dof,1); % acceleration SDB
% components at each gauss point
ng = 2;
pg = ng^2;
Eps=zeros(3*pg,size(TB,1));
eps=zeros(3*pg,size(TB,1));
sig=zeros(3*pg,size(TB,1));
Sig=zeros(3*pg,size(TB,1));
P_past = zeros(nnB*dof,1);
%% Calculate the matrices for the calculation of the internal forces
[Mat_epsi_n,Mat_Sig_n,vect_epsi_n,vect_epsi_n_2,vect_sig_n,vect_Epsi_n,...
    vect_Sig_n,vect_past,Mat_Bepsi_n,Mat_BQ_n,Mat_Fepsi_n,Mat_FQ_n,...
    Mat_Dsig,Mat_Epsi_n] = calculate_cste_matrices(Nef_2, ng,k_1,mu_1,TB,XB,...
    a0,b0,L2,n,x0,h,cs,b_scal,dt2);

%% Plot
% Reshape for plotting
subplot(2,2,1);
coord1=0:Le_1:L1;
coord2=L1:Le_2:L2+L1;
axis([0 (L1+L2) -5 5]);
set(gca,'NextPlot','replacechildren');
size1= size(coord1,2);
plot(coord1,u1(1:2:102),'k-',coord2,u2(1:2:82),'r-','linewidth',2);

subplot(2,2,2);
coord1=0:Le_1:L1;
coord2=L1:Le_2:L2+L1;
axis([0 (L1+L2) -5 5]);
set(gca,'NextPlot','replacechildren');
plot(coord1,u1(2:2:102),'k-',coord2,u2(2:2:82),'r-','linewidth',2);

subplot(2,2,[3,4]);
% init
XA_new = zeros(nnA/2+nnB/2,2);
X_plot = zeros(nnA/2+nnB/2,1); Y_plot = zeros(nnA/2+nnB/2,1);
% SDA
X_plot(1:51,1) = XA(1:51,1); X_plot(1:51,2) = XA(52:102,1);
Y_plot(1:51,1) = XA(1:51,2); Y_plot(1:51,2) = XA(52:102,2);
% SDB
X_plot(52:92,1) = XB(1:41,1); X_plot(52:92,2) = XB(42:82,1);
Y_plot(52:92,1) = XB(1:41,2); Y_plot(52:92,2) = XB(42:82,2);
surf(X_plot,Y_plot,ones(size(X_plot)));
view(2);
axis([0 (L1+L2) 0 h -1 5]);
axis(gca,'on');
caxis([-1 1]);
colorbar
F(ntps_1) = struct('cdata',[],'colormap',[]);
F(1) = getframe(gcf);

t1=0;t2=0;
Forn1 = zeros(dof*nnA,1);
Forn1(1,1) = 1;  Forn1(103,1) = 1; 

for tt = 2:ntps_1
    % Calculate internal forces
    [vect_past, P_past] = internal_forces_PML(dof,ng,Mat_Sig_n,vect_Sig_n,Mat_epsi_n,...
        vect_epsi_n,Mat_Epsi_n,vect_Epsi_n,P_past,TB,Nef_2);
    
    % Ricker force
    f1 = F_ricker(1E6,t1); % force
    Fext1 = f1*Forn1;      % application
    % Predictors
    t1=t1+dt1;
    up1 = u1+(dt1*v1)+(((0.5-beta1)*(dt1^2))*a1);  
    vp1 = v1+(((1.0-gamma1)*dt1)*a1);
    % Righthand term
    rhs = Fext1-(KA*up1); 
    % acceleration
    a1 = rhs./diagM;
    % update
    v1  = vp1+((gamma1*dt1)*a1); 
    u1  = up1+((beta1*(dt1^2))*a1); 
    
    v2(1) = v1(101);
    v2(83) = v1(203);
    a2(1) = a1(101);
    a2(83) = a1(203);
    
    % Predictors on SDB
    t2=t2+dt2; % increment
    up2 = u2 + dt2*v2 + (dt2)^2*(0.5-beta2)*a2;
    
    vp2 = v2+dt2*(1-gamma2)*a2;
    % calculation of acceleration
    P_past= P_past + (CCB+CB)*vp2 + (KKB+KB)*up2 ;
    a2= (inv_M_eff) * (-P_past);
    % update
    u2 = beta2*(dt2)^2*a2+up2;
    v2 = gamma2*dt2*a2+vp2;
  
    v1(101) = v2(1);
    v1(203) = v2(83);
    a1(101) = a2(1);
    a1(203) = a2(83);

    
    % Junction
    % Continuity of displacement, velocity and acceleration
    
    % Update of physical quantities
    [P_past,vect_epsi_n_2,vect_epsi_n,vect_sig_n,vect_Epsi_n,...
    vect_Sig_n] = update_quantities(dt2,nnB,dof,ng,Nef_2,global_pos_elemB,u2...
    ,v2,vect_epsi_n_2,vect_epsi_n,Mat_Bepsi_n,Mat_BQ_n,Mat_Fepsi_n,...
    Mat_FQ_n,vect_sig_n,Mat_Dsig,vect_Epsi_n,vect_Sig_n);
    
    
    
    
    %% Plot
    subplot(2,2,1);
    plot(coord1,u1(1:2:102),'k-',coord2,u2(1:2:82),'r-','linewidth',2);
    subplot(2,2,2);
    plot(coord1,u1(2:2:102),'k-',coord2,u2(2:2:82),'r-','linewidth',2);
%     subplot(2,2,[3,4]);
%     XA_new(1:102,1) = XA(:,1) + u1(1:2:204);
%     XA_new(1:102,2) = XA(:,2) + u1(2:2:204);
%     XA_new(103:184,1) = XB(:,1) + u2(1:2:164);
%     XA_new(103:184,2) = XB(:,2) + u2(2:2:164);
%     Z(1:51,1) = u1(1:2:102); Z(1:51,2) = u1(103:2:204);
%     Z(52:92,1) = u2(1:2:82); Z(52:92,2) = u2(83:2:164);
%     % SDA
%     X_plot(1:51,1) = XA_new(1:51,1); X_plot(1:51,2) = XA_new(52:102,1);
%     Y_plot(1:51,1) = XA_new(1:51,2); Y_plot(1:51,2) = XA_new(52:102,2);
%     % SDB
%     X_plot(52:92,1) = XA_new(103:143,1); X_plot(52:92,2) = XA_new(144:184,1);
%     Y_plot(52:92,1) = XA_new(103:143,2); Y_plot(52:92,2) = XA_new(144:184,2);
% 
%     surf(X_plot,Y_plot,Z);
%     axis([0 (L1+L2) 0 h -1 5]);
%     caxis([-1 1]);
%     colorbar
%     shading interp
%     view(2)
    F(tt)=getframe;
    % time step
    tt
end







