%% 1 dimentional SDOF bar element of PML Stability

clear all; close all;
% Newmark parameters
gamma = 1/2; beta = 1/4; % Implicit
%gamma = 1/2; beta = 0;   % Explicit

% parameters of PML
rho = 1700; % mass density
E = 10e5;   % Young modulus
nu = 0.24; % Poisson modulus
mu = E/(2*(1+nu)); % 1st Lame Coefficient
cs = (mu/rho)^0.5 ; % Shear wave velocity
S = 1;      % area of the bar
L = 200;    % length of PML
x0 = 250;   % start of the PML: end of phsical medium
Le = 1;     % length of the element
% shape function and its derivative
N1 = @(xi) 1/2*(1-xi);
N2 = @(xi) 1/2*(1+xi);
N = @(xi)[N1(xi) N2(xi)];
R=(1/Le)*[-1. 1.];
% attenuation functions
Rpp = 10^-3;  % reflection wanted
n = 2;        % order of attenuation
alpha_kucu = sqrt(E/rho)*(n+1)/(2*L)*log(1/Rpp); % attenuation ...
%                                  coefficient calculated by Kucukcoban
% alpha_kucu = 0; % no damping
a0 = alpha_kucu;
b0 = alpha_kucu;
fe = @(x) a0*((x-x0)/L)^n; 
fp = @(x) b0*((x-x0)/L)^n; 
nbx_pml =1; % number of element
% gauss quadrature: 2 gauss points
xi = [-1/sqrt(3) 1/sqrt(3)];
w = [1 1];
J=Le/2;

% Calculation of omega
M = rho*S*( (w(1)*N(xi(1))'*N(xi(1))*J)+ (w(2)*N(xi(2))'*N(xi(2))*J) ); % mass matrix
K = E*(w(1)*R'*R*J + w(2)*R'*R*J);
% K(1,1) = 0; K(1,2) = 0; K(2,1) = 0; % encastred beam  
eigval = eig(inv(M)*K);
omega= sqrt(eigval(2));
h=[0.001:0.001:10];
Omega = sqrt(eigval(2))*h;
l_Omega = Omega - 4.0 ;
[val_min,i_min] = min(abs(l_Omega)) ;

me = zeros(2,2); % Submatrix
M = zeros(2,2);
for e=1:nbx_pml
    % Position of x in physical space
    e1 = Le*(e-1) + L;
    e2 = Le*(e) + L;
    % Position of gauss point in physical space
    pos1 = 1/2*(1-xi(1))*e1 + 1/2*(1+xi(1))*e2;
    pos2 = 1/2*(1-xi(2))*e1 + 1/2*(1+xi(2))*e2;
    
    me = rho*Le/2*S*( (w(1)*N(xi(1))'*N(xi(1)))+ (w(2)*N(xi(2))'*N(xi(2))) );
    M(e:e+1,e:e+1) = M(e:e+1,e:e+1) + me(:,:);
end

ce = zeros(2,2);
C = zeros(2,2);
for e=1:nbx_pml
    % Position of x in physical space
    e1 = Le*(e-1) + x0;
    e2 = Le*(e) + x0;
    % Position of gauss point in physical space
    pos1 = 1/2*(1-xi(1))*e1 + 1/2*(1+xi(1))*e2;
    pos2 = 1/2*(1-xi(2))*e1 + 1/2*(1+xi(2))*e2;
    ce =    (w(1) * (fp(pos1))*[N1(xi(1))^2 N1(xi(1))*N2(xi(1)) ; N1(xi(1))*N2(xi(1)) N2(xi(1))^2] ...
        + w(2) * (fp(pos2))*[N1(xi(2))^2 N1(xi(2))*N2(xi(2)); N1(xi(2))*N2(xi(2)) N2(xi(2))^2]);
    C(e:e+1,e:e+1) = C(e:e+1,e:e+1) + rho * cs *S* Le/2 * ce(:,:); % Pourquoi pas divis� par la longueur de la pml
end

for i=1:length(h)
    dt = h(i);
    ke = zeros(2,2);
    K = zeros(2,2);
    for e=1:nbx_pml
        % Position of x in physical space
        e1 = Le*(e-1) + x0;
        e2 = Le*(e) + x0;
        % Position of gauss point in physical space
        pos1 = 1/2*(1-xi(1))*e1 + 1/2*(1+xi(1))*e2;
        pos2 = 1/2*(1-xi(2))*e1 + 1/2*(1+xi(2))*e2;
        % Calculation oe 1/alpha
        alpha_inv = w(1)/((1+fe(pos1))+cs/L*fp(pos1)*dt) + w(2)/((1+fe(pos2))+cs/L*fp(pos2)*dt); 
        ke = E * alpha_inv * (R)'*R*S* Le/2;
        % Assembling
        K(e:e+1,e:e+1) = K(e:e+1,e:e+1) + ke(:,:); 
        Ap(1,1) = 1/2*cs/L*E*S*(w(1)*fp(pos1)/((1+fe(pos1))+cs/L*fp(pos1)*dt))*(-1);
        Ap(1,2) = 1/2*cs/L*E*S*(w(1)*fp(pos2)/((1+fe(pos2))+cs/L*fp(pos2)*dt))*(-1);
        Ap(2,1) = 1/2*cs/L*E*S*(w(1)*fp(pos1)/((1+fe(pos1))+cs/L*fp(pos1)*dt));
        Ap(2,2) = 1/2*cs/L*E*S*(w(1)*fp(pos2)/((1+fe(pos2))+cs/L*fp(pos2)*dt));
    end
    tempB = 1/Le.*[-1 1; -1 1];
    A = [eye(size(M^-1*K))+h(i)*(1-gamma)*M^-1*C gamma*h(i)*M^-1*K -h(i)*gamma*M^-1*Ap;
         h(i)^2*beta*M^-1*C eye(size(M^-1*K))+h(i)^2*beta*M^-1*K -h(i)^2*beta*M^-1*Ap;
         zeros(size(M)) -h(i)*tempB*alpha_inv eye(2)];
     B = -[-eye(size(M^-1*K))+ h(i)*(1-gamma)*M^-1*C (1-gamma)*h(i)*M^-1*K -h(i)*(1-gamma)*M^-1*Ap;
         -eye(size(M^-1*K))*h(i)+h(i)^2*(1/2-beta)*M^-1*C -eye(size(M^-1*K))+h(i)^2*(1/2-beta)*M^-1*K -h(i)^2*(1/2-beta)*M^-1*Ap;
         zeros(2) zeros(2) -eye(2)+h(i)*cs/L*( w(1)*fp(pos1)/((1+fe(pos1))+cs/L*fp(pos1)*dt) + w(2)*fp(pos2)/((1+fe(pos2))+cs/L*fp(pos2)*dt))];
    % Calculation of the amplification matrix
    A_xi(:,:,i) = A^-1 * B;
end
[V,D]=eigenshuffle(A_xi);
for i = 1:length(h)
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
for i = 1:length(h)
    if mod(i,20)==0
%         plot(eigen_val((i-1)*2+1,:)...
%             ,eigen_val((i-1)*2+2,:),'+');
        plot(eigen_val((i-1)*2+1,5)...
            ,eigen_val((i-1)*2+2,5),'+r');
        plot(eigen_val((i-1)*2+1,6)...
            ,eigen_val((i-1)*2+2,6),'og');
        hold on ;
    end
end
grid on ;
axis equal ;
axis([-2 2 -2 2]) ;

% spectral radius
figure(2) ;
clf ;
for i = 1:length(h)
    if mod(i,20)==0 
        %semilogx(Omega(i),lambda_xi(i,:),'+') ;
        semilogx(Omega(i),lambda_xi(i,5),'+r') ;
        semilogx(Omega(i),lambda_xi(i,6),'og') ;
        hold on ;
    end
end
grid on ;
% 
% Periodicity error
figure(3) ;
clf ;
for i = 1:i_min
    % plot(Omega(i),PE_xi(i,:),'+') ;
    plot(Omega(i),PE_xi(i,5),'+r') ;
    plot(Omega(i),PE_xi(i,6),'og') ;
    hold on ;
end
grid on ;

% numerical damping ratio
figure(4) ;
clf ;
for i = 1:i_min
    % plot(Omega(i),AD_xi(i,:),'+') ;
    plot(Omega(i),AD_xi(i,5),'+r') ;
    plot(Omega(i),AD_xi(i,6),'og') ;
    hold on ;
end
grid on ;
% 
