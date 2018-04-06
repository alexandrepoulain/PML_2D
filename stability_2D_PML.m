%% Numerical stability of the 2D PML
% Stability of the integration method
% Author: Alexandre Poulain
% Supervisor: Michael Brun
% On 1 linear quadrilateral (4-noded) element
clear all; close all;

% Newmark parameters
gamma = 0.5; beta = 0.25;       % Implicit
% gamma = 0.5; beta = 0;        % Explicit 

% Time stepping parameters
omega=2*pi;                     % pulsation
dt=[0.001:0.1:10];            % time step
Omega = omega*dt ;              % frequency

% Parameters of the material
rho = 1700;                     % density
nu = 0.24;                      % Poisson modulus
E = 1e07;                       % Young modulus
mu = E/(2*(1+nu));              % 1st Lame Coefficient
lmbd = E*nu/((1+nu)*(1-2*nu));  % 2nd Lame Coefficient
bulk = E/3/(1-2*nu);            % Bulk modulus
cs = (mu/rho)^0.5 ;             % Shear wave velocity
cp = sqrt((lmbd+2*mu)/rho);     % Pressure wave velocity
b_scal = cs;                    % scalar
L1 = 250;                       % length of medium A

% Assembly element matrices for the 2D PML: M, C, Ceff, K, Keff
Nex = 1; Ney = 1;             % One element in x and y direction
Lex = 5; Ley = 5;             % Length of an element in x and y direction
[XB,TB] = ...
    topology(0,0,Lex,Ley,Nex,Ney); % Coordinates and connectivity matrix 
                       %(The element of PML starts at x_start = y_start= 0)
nnodes = 4;                     % number of nodes in an element
dof = 2;                        % number of degrees of freedom
ng = 2;                         % order of quadrature

% Parameters for the PML
L2 = 200; h = 5;                % length of PML, length in y direction 
n = 2;                          % order for attenuation functions
Rpp = 1e-3;                     % reflexion wanted
alpha_kucu = sqrt(E/rho)...
    *(n+1)/(2*L2)*log(1/Rpp);   % attenuation coefficient (ref:Kucukcoban)
a0 = [alpha_kucu,alpha_kucu];   % coefficient attenuation evanescent waves
b0 = [alpha_kucu,alpha_kucu];   % coefficient attenuation propagating waves
x0 = L1;                       % start of the PML
% From local to global and in 2D
global_pos_elemB = zeros(ng^2*dof,size(TB,1));
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
MB = assembleM_PML(rho,XB,TB,a0,x0,L2,h,n,dof);
% stiffness
KB = assembleK_PML(rho,cs,b_scal,XB,TB,a0,b0,x0,L2,h,n,dof);
% damping
CB = assembleC_PML(rho,cs,b_scal,XB,TB,a0,b0,x0,L2,h,n,dof);
for i=1:length(dt)
    i
    % Calculation of these matrices here because they depend on dt 
    % effective stiffness
    KKB = assemble_effKB_PML(bulk,mu,cs,b_scal,dt(i),XB,TB,a0,b0,x0,L2,h,n,dof);
    % effective damping
    CCB = assemble_effCB_PML(bulk,mu,cs,b_scal,dt(i),XB,TB,a0,b0,x0,L2,h,n,dof);
    
    % Calculation of Amplification matrices
    A = calculation_amplification_A(dt(i),beta,gamma,bulk,mu,ng,TB,...
    XB,L2,n,a0,b0,x0,h,cs,b_scal,MB,CB,CCB,KB,KKB,nnodes,dof); 
    B = calculation_amplification_B(dt(i),beta,gamma,bulk,mu,ng,TB,...
    XB,L2,n,a0,b0,x0,h,cs,b_scal,MB,CB,CCB,KB,KKB,nnodes,dof);
    % Calculation of the amplification matrix
    A_xi = A^-1 * B;
    % calculation of its eigenvalues (number of eig_val depends on ng)
    d_xi(i,:)=eig(A_xi);
    % For each eigenvalue
    for k = 1:size(d_xi(i,:),2)
       eigen_val((i-1)*2+1,k) =  real(d_xi(i,k)); % real part
       eigen_val((i-1)*2+2,k) =  imag(d_xi(i,k)); % imaginary part
       % spectral radius
       lambda_xi(i,k) = abs(d_xi(i,k));
       % relative periodicity error
       PE_xi(i,k) = ( omega * dt(i) / abs(angle(d_xi(i,k))) ) - 1;
       % algorithmic damping ratio
       AD_xi(i,k)= -log(lambda_xi(i,k)) / abs(angle(d_xi(i,k)));
    end
end
% plots
% eigen values in imaginary plan
figure(1) ;
clf ;
for i = 1:length(dt)
    for k = 1:length(d_xi(i,:)) 
        if mod(i,20)==0 
            plot(eigen_val((i-1)*2+1,:)...
                ,eigen_val((i-1)*2+2,:),'+');
            hold on ;
        end
    end
end
grid on ;
axis equal ;
axis([-2 2 -2 2]) ;






