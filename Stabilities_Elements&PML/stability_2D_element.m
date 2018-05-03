%% Stability of Newmark integrator 
% Stability of the integration method
% Author: Alexandre Poulain
% Supervisor: Michael Brun
% On 1 linear quadrilateral (4-noded) element
clear all; close all;

% Newmark parameters
gamma = 0.5; beta = 0.25;       % Implicit
gamma = 0.5; beta = 0;        % Explicit 

% Parameters of the material
rho = 1700;                     % density
nu = 0.24;                      % Poisson modulus
E = 1e07;                       % Young modulus
mu = E/(2*(1+nu));              % 1st Lame Coefficient
lmbd = E*nu/((1+nu)*(1-2*nu));  % 2nd Lame Coefficient
bulk = E/3/(1-2*nu);            % Bulk modulus
cs = (mu/rho)^0.5 ;             % Shear wave velocity
cp = sqrt((lmbd+2*mu)/rho);     % Pressure wave velocity

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
% lumped Mass
NA = dof * nnodes ;
diagM = zeros(NA,1) ;
for i=1:NA 
  diagM(i) = sum(M(i,1:NA)) ;
end
% Only if explicit 
M = diag(diagM) ;
% stiffness
K = assembleK(XA,TA,dof,E,nu,Nef);

% Modes 
[V0,D0] = eigenshuffle(K);
% Calculation of iniial vibration modes
[init_vec,init_modes] = eigenshuffle(M^-1*K);

% Time stepping parameters
omega= max(sqrt(init_modes));                    % pulsation
dt=[[1e-4:1e-4:1e-1]];            % time step
Omega = omega*dt ;             % frequency
l_Omega = Omega - 8.0 ;
[val_min,i_min] = min(abs(l_Omega)) ;
for i=1:length(dt)
    A = [eye(size(M^-1*K)) gamma*dt(i)*M^-1*K;
         zeros(size(M^-1*K)) eye(size(M^-1*K))+dt(i)^2*beta*M^-1*K];
     B = -[-eye(size(M^-1*K)) (1-gamma)*dt(i)*M^-1*K;
         -eye(size(M^-1*K))*dt(i) -eye(size(M^-1*K))+dt(i)^2*(1/2-beta)*M^-1*K];
    % Calculation of the amplification matrix
    A_xi(:,:,i) = inv(A) * B;
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
       PE_xi(i,k) = ( Omega(i) / abs(angle(D(k,i))) ) - 1;
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
%         plot(eigen_val((i-1)*2+1,3)...
%             ,eigen_val((i-1)*2+2,3),'r+');
%          plot(eigen_val((i-1)*2+1,4)...
%              ,eigen_val((i-1)*2+2,4),'og');
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
%         semilogx(Omega(i),lambda_xi(i,12),'+b') ;
%         semilogx(Omega(i),lambda_xi(i,14),'+b') ;
%         semilogx(Omega(i),lambda_xi(i,15),'+b') ;
        
%         semilogx(Omega(i),lambda_xi(i,3),'+r') ;
%         semilogx(Omega(i),lambda_xi(i,4),'og') ;
%         hold on ;
%    
% end
semilogx(Omega(1:1:length(Omega)),lambda_xi(1:1:length(Omega),:),'.k') ;

title('spectral radius')

grid on ;

%% Periodicity error
figure(3) ;
clf ;
for i = 1:i_min
    
    plot(Omega(i),PE_xi(i,12),'+r') ;
    plot(Omega(i),PE_xi(i,14),'+y') ;
    plot(Omega(i),PE_xi(i,15),'+b') ;
%     plot(Omega(i),PE_xi(i,3),'+r') ;
%     plot(Omega(i),PE_xi(i,4),'og') ;
    hold on ;
end
title('relative periodicity error')

grid on ;


% numerical damping ratio
figure(4) ;
clf ;
for i = 1:i_min
    plot(Omega(i),AD_xi(i,12),'+r') ;
    plot(Omega(i),AD_xi(i,14),'+y') ;
    plot(Omega(i),AD_xi(i,15),'+b') ;
%     plot(Omega(i),AD_xi(i,3),'+r') ;
%     plot(Omega(i),AD_xi(i,4),'og') ;
    hold on ;
end
title('Numerical damping ratio')

grid on ;
