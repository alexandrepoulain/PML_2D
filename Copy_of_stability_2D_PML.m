%% Numerical stability of the 2D PML
% Stability of the integration method
% Author: Alexandre Poulain
% Supervisor: Michael Brun
% On 1 linear quadrilateral (4-noded) element
clear all; close all;

% Newmark parameters
gamma = 0.5; beta = 0.25;       % Implicit
gamma = 0.5; beta = 0;        % Explicit 
%x domaine
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
L1 = 1;                       % length of medium A

% Assembly element matrices for the 2D PML: M, C, Ceff, K, Keff
Nef = 1;                      % number of elements
Nex = 1; Ney = 1;             % One element in x and y direction
Lex = 1; Ley = 1;             % Length of an element in x and y direction
[XA,TA] = ...
    topology(0,0,Lex,Ley,Nex,Ney); % Coordinates and connectivity matrix 
                       %(The element of PML starts at x_start = y_start= 0)
nnodes = 4;                     % number of nodes in an element
dof = 2;                        % number of degrees of freedom
ng = 2;                         % order of quadrature

% assemble mass and stiffness matrices
% mass 
M = assembleM(XA,TA,dof,Nef,rho);
% % lumped Mass
% NA = dof * nnodes ;
% diagM = zeros(NA,1) ;
% for i=1:NA 
%   diagM(i) = sum(M(i,1:NA)) ;
% end
% %Only if explicit 
% M = diag(diagM) ;
% stiffness
K = assembleK(XA,TA,dof,E,nu,Nef);

% Calculation of iniial vibration modes
[init_vec,init_modes] = eigenshuffle(M^-1*K);

% Time stepping parameters
omega= max(sqrt(init_modes));                    % pulsation

% Time stepping parameters
dt=[[1e-3:1e-3:1]];            % time step
Omega = omega*dt ;              % frequency
l_Omega = Omega - 8.0 ;
[val_min,i_min] = min(abs(l_Omega)) ;



% Assembly element matrices for the 2D PML: M, C, Ceff, K, Keff
Nex = 1; Ney = 1;             % One element in x and y direction
Lex = 1; Ley = 1;             % Length of an element in x and y direction
[XB,TB] = ...
    topology(Lex,0,2*Lex,Ley,Nex,Ney); % Coordinates and connectivity matrix 
                       %(The element of PML starts at x_start = y_start= 0)
nnodes = 4;                     % number of nodes in an element
dof = 2;                        % number of degrees of freedom
ng = 2;                         % order of quadrature

% Parameters for the PML
L2 = 1; h = 1;                % length of PML, length in y direction 
n = 2;                          % order for attenuation functions
Rpp = 1e-1;                     % reflexion wanted
alpha_kucu = sqrt(E/rho)...
    *(n+1)/(2*L2)*log(1/Rpp);   % attenuation coefficient (ref:Kucukcoban)
a0 = [alpha_kucu,alpha_kucu];   % coefficient attenuation evanescent waves
b0 = [alpha_kucu,alpha_kucu];   % coefficient attenuation propagating waves
a0 = [0,0];
b0 = [0,0]; 
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
% lumped Mass
% NA = dof * nnodes ;
% diagM = zeros(NA,1) ;
% for i=1:NA 
%   diagM(i) = sum(MB(i,1:NA)) ;
% end
% % Only if explicit 
% MB = diag(diagM) ;
% stiffness
KB = assembleK_PML(rho,cs,b_scal,XB,TB,a0,b0,x0,L2,h,n,dof);
% damping
CB = assembleC_PML(rho,cs,b_scal,XB,TB,a0,b0,x0,L2,h,n,dof);
%% Calculation
for i=1:length(dt)
    i
    % Calculation of these matrices here because they depend on dt 
    % effective stiffness
    KKB = assemble_effKB_PML(bulk,mu,cs,b_scal,dt(i),XB,TB,a0,b0,x0,L2,h,n,dof);
    % effective damping
    CCB = assemble_effCB_PML(bulk,mu,cs,b_scal,dt(i),XB,TB,a0,b0,x0,L2,h,n,dof);
    
    % Calculation of Amplification matrices
    A = Copy_of_calculation_amplification_A(dt(i),beta,gamma,bulk,mu,ng,TB,...
    XB,L2,n,a0,b0,x0,h,cs,b_scal,MB,CB,CCB,KB,KKB,nnodes,dof); 
    B = Copy_of_calculation_amplification_B(dt(i),beta,gamma,bulk,mu,ng,TB,...
    XB,L2,n,a0,b0,x0,h,cs,b_scal,MB,CB,CCB,KB,KKB,nnodes,dof);
    % Calculation of the amplification matrix
    A_xi(:,:,i) = A^(-1) * (-B);

end
[V,D]=eigenshuffle(A_xi);

for i = 1:length(dt)
       % calculation of its eigenvalues (number of eig_val depends on ng)
    %d_xi(i,:)=eig(A_xi(:,:,i));
    % For each eigenvalue
    for k = 1:size(D(:,i),1)
       eigen_val((i-1)*2+1,k) =  real(D(k,i)); % real part
       eigen_val((i-1)*2+2,k) =  imag(D(k,i)); % imaginary part
       % spectral radius
       lambda_xi(i,k) = abs(D(k,i));
       % relative periodicity error
       PE_xi(i,k) = ( omega * dt(i) / abs(angle(D(k,i))) ) - 1;
       % algorithmic damping ratio
       AD_xi(i,k)= -log(lambda_xi(i,k)) / abs(angle(D(k,i)));
    end 
end
%% plots
% eigen values in imaginary plan
figure(1) ;
clf ;
for i = 1:length(dt)
    if mod(i,20)==0 
        plot(eigen_val((i-1)*2+1,:)...
            ,eigen_val((i-1)*2+2,:),'+');
        hold on ;
    end
end
grid on ;
title('Eigenvalues in imaginary plan')

axis equal ;
axis([-2 2 -2 2]) ;

% spectral radius
figure(2) ;
clf ;
% for i = 1:length(dt)
%     if mod(i,20)==0 
        %semilogx(Omega(i),lambda_xi(i,:),'+') ;
        
%         hold on ;
%     end
% end
semilogx(Omega(1:1:length(Omega)),lambda_xi(1:1:length(Omega),:),'.k') ;
title('spectral radius')

grid on ;

%% Periodicity error
tic
figure(3) ;
clf ;
for i = 1:i_min
    plot(Omega(i),PE_xi(i,25:26),'+') ;
    plot(Omega(i),PE_xi(i,21:22),'+') ;
    hold on ;
end
title('relative periodicity error')

grid on ;


% numerical damping ratio
figure(4) ;
clf ;
for i = 1:i_min
    plot(Omega(i),AD_xi(i,25:26),'+') ;
    plot(Omega(i),AD_xi(i,21:22),'+') ;
    hold on ;
end
title('Numerical damping ratio')

grid on ;

toc

