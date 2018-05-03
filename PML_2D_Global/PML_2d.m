clear; clc; close all;
%% PML 2D enhanced (2nd order) bar Test
% Parameters of input Ricker Wave
tp = 3; % Fundamental period
ts = 3; % Time shift
A = 1; % Amplitude


%% Material parameters
%%%%%%%%%%%%%%% SDA %%%%%%%%%%%%%%%
rho = 1700; % density
nu = 0.24; % Poisson modulus
E = 1e07; % Young modulus
mu = E/(2*(1+nu)); % 1st Lam� Coefficient
lmbd = E*nu/((1+nu)*(1-2*nu)); % 2nd Lam� Coefficient
bulk = E/3/(1-2*nu);   % Bulk modulus
cs = (mu/rho)^0.5 ; % Shear wave velocity
cp = sqrt((lmbd+2*mu)/rho); % Pressure wave velocity
lmbd_s = cs*tp; % wave length of shear wave
b_scal = cs;
L1 = 250; % Size of the soil
h = 5; % width of the plane
Lef = 5;
Nef = floor(L1 / Lef) ; % number of finite element in SDA
Nex1 = Nef; % nb elem in x direction
Ney1 = 1; % nb elem in y direction
% Space discretization (4-nodes bilinear elements)
Le_1 = 5; % Length of element 
dof = 2; % Per element
% Square elements
a = Lef;
b = Lef;
%%%%%%%%%%%%%%% SDB %%%%%%%%%%%%%%%
L2 = 200;   % Length of the PML
Le_2 = 5;   % Length of an element
Nef_2 = floor( L2 / Le_2 ); % Number of elements
Nef_tot = Nef + Nef_2 ;
L_tot = L1+L2;
% Number of nodes
Nex2 = Nef_2 ; % number of nodes in x-direction
Ney2 = 1 ; % number of nodes in y-direction
n = 2; % order for attenuation functions
Rpp = 1e-3; % reflexion wanted
alpha_kucu = sqrt(E/rho)*(n+1)/(2*L2)*log(1/Rpp); % attenuation coefficient calculated by Kucukcoban
% alpha_kucu = 0;
% The coefficient along  y-direction is null 
a0 = [alpha_kucu,0]; % coefficient of attenuation evanescent waves
b0 = [alpha_kucu,0]; % coefficient of attenuation propagating waves
x0 = L1; % start of the PML




%% global Topology
% construction of the parameters arrays
param_med = {'MED', rho, nu, E, mu, lmbd, bulk, cs, cp, lmbd_s,...
    b_scal, L1, h, Nef, Nex1, Ney1, Le_1, dof};
param_PML = {'PML', rho, nu, E, mu, lmbd, bulk, cs, cp, lmbd_s,...
    b_scal, L2, h, Nef_2, Nex2, Ney2, Le_2, dof, a0, b0, x0};

[X,T,G] = global_topology(L1,h,Nex1,Ney1,Le_1,L2,Nex2,Ney2,Le_2,param_med,param_PML);

% From local to global and in 2D
global_pos_elemB = zeros(4*dof,size(T,1));
for i=1:size(T,1)
    nnodes = size(T,2)-1;  
    for s = 1:nnodes
        for j = 1:dof
          global_pos_elemB((s-1)*dof+j,i) = (T(i,s)-1)*dof + j ;
        end
    end
end
%% Time evolution
% parameters
dt_CFL = min(a,b) / sqrt(E/rho) ;
dt = dt_CFL*0.8;
t_fin = 100.;
ntps_1 = round(t_fin/dt);
gamma=0.5; beta=0.25;    % Implicit Newmark
% gamma=0.5; beta=0.0;    % Explicit Newmark

%% Construct global matrices
ng=2;
[Mass, Stiff, Damp] = construct_global(X,T,G,rho,cs,b_scal,bulk,mu,L2,x0,h,n,dt,dof,ng,a0,b0);
% NN = dof * size(X,1) ;
% diagM = zeros(NN,1) ;
% for i=1:NN 
%   diagM(i) = sum(Mass(i,1:NN)) ;
% end
Mass_eff = sparse(Mass+gamma*dt*Damp+ beta*(dt)^2*Stiff);
inv_M_eff=inv(Mass_eff);

%% Initializing
nn = size(X,1);         % number of nodes 
u = zeros(nn*dof,1); % displacement 
v = zeros(nn*dof,1);       % velocity 
a = zeros(nn*dof,1);       % acceleration 

% components at each gauss point
ng = 2;                         % number of gauss point in each directions
pg = ng^2;                      % number of gauss points in total
Eps=zeros(3*pg,size(T,1));  % Integral of strain
eps=zeros(3*pg,size(T,1));  % strain
sig=zeros(3*pg,size(T,1));  % stress
Sig=zeros(3*pg,size(T,1));  % integral of stress
P_past = zeros(nn*dof,1);       % Vector of Internal forces
%% Calculate the matrices for the calculation of the internal forces

[Mat_epsi_n,Mat_Sig_n,vect_epsi_n,vect_epsi_n_2,vect_sig_n,vect_Epsi_n,...
    vect_Sig_n,vect_past,Mat_Bepsi_n,Mat_BQ_n,Mat_Fepsi_n,Mat_FQ_n,...
    Mat_Dsig,Mat_Epsi_n] = calculate_cste_matrices(Nef_tot, ng,bulk,mu,T,X,G,...
    a0,b0,L2,n,x0,h,cs,b_scal,dt);

%% Plot
% Reshape for plotting
subplot(2,2,1);
coord1=0:Le_1:L1;
coord2=L1:Le_2:L2+L1;
axis([0 (L1+L2) -5 5]);
set(gca,'NextPlot','replacechildren');
size1= size(coord1,2);
plot(coord1,u(1:2:102),'k-',coord2,u(102:2:182),'r-','linewidth',2);

subplot(2,2,2);
coord1=0:Le_1:L1;
coord2=L1:Le_2:L2+L1;
axis([0 (L1+L2) -5 5]);
set(gca,'NextPlot','replacechildren');
plot(coord1,u(2:2:102),'k-',coord2,u(101:2:182),'r-','linewidth',2);

% subplot(2,2,[3,4]);
% % init
% XA_new = zeros(nnA/2+nnB/2,2);
% X_plot = zeros(nnA/2+nnB/2,1); Y_plot = zeros(nnA/2+nnB/2,1);
% % SDA
% X_plot(1:92,1) = X(1:92,1); X_plot(1:92,2) = XA(52:102,1);
% Y_plot(1:92,1) = X(1:51,2); Y_plot(1:92,2) = XA(52:102,2);
% % SDB
% X_plot(52:92,1) = XB(1:41,1); X_plot(52:92,2) = XB(42:82,1);
% Y_plot(52:92,1) = XB(1:41,2); Y_plot(52:92,2) = XB(42:82,2);
% surf(X_plot,Y_plot,ones(size(X_plot)));
% view(2);
% axis([0 (L1+L2) 0 h -1 5]);
% axis(gca,'on');
% caxis([-1 1]);
% colorbar
% F(ntps_1) = struct('cdata',[],'colormap',[]);
% F(1) = getframe(gcf);

t=0;
Forn = zeros(dof*nn,1);
Forn(1,1) = 1;  Forn(183,1) = 1; 

for tt = 2:ntps_1
    % Calculate internal forces
    [vect_past, P_past] = internal_forces_PML(dof,ng,Mat_Sig_n,vect_Sig_n,Mat_epsi_n,...
        vect_epsi_n,Mat_Epsi_n,vect_Epsi_n,P_past,T,G,Nef_tot);
    
    % Ricker force
    f1 = F_ricker(1E6,t); % force
    Fext = f1*Forn;      % application
    % Predictors
    t=t+dt;
    up = u+(dt*v)+(((0.5-beta)*(dt^2))*a);  
    vp = v+(((1.0-gamma)*dt)*a);
    % acceleration
%    P_past = P_past + Stiff*vp + Damp*up ;
%    a=   (Mass^(-1))*(Fext-P_past);
    rhs = Fext-(P_past+Damp*vp+Stiff*up); 
    % acceleration
    a = (inv_M_eff) * rhs;
    % update
    v  = vp+((gamma*dt)*a); 
    u  = up+((beta*(dt^2))*a); 
    if tt >= 80 
        "here"
    end
    % Update of physical quantities
    [P_past,vect_epsi_n_2,vect_epsi_n,vect_sig_n,vect_Epsi_n,...
    vect_Sig_n] = update_quantities(G,dt,nn,dof,ng,Nef_tot,global_pos_elemB,u...
    ,v,vect_epsi_n_2,vect_epsi_n,Mat_Bepsi_n,Mat_BQ_n,Mat_Fepsi_n,...
    Mat_FQ_n,vect_sig_n,Mat_Dsig,vect_Epsi_n,vect_Sig_n);
    
    
    
    
    %% Plot
    subplot(2,2,1);
    plot(coord1,u(1:2:102),'k-',coord2,u(102:2:182),'r-','linewidth',2);
    subplot(2,2,2);
    plot(coord1,u(2:2:102),'k-',coord2,u(101:2:182),'r-','linewidth',2);
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







