%% 1 dimentional SDOF bar element Stability

clear all; close all;
% Newmark parameters
gamma = 1/2; beta = 1/4; % Implicit
gamma = 1/2; beta = 0;   % Explicit
rho = 1700; % mass density
E = 10e3;
A = 1; % area of the bar
Le = 1; % length of the element
% shape function and its derivative
N1 = @(xi) 1/2*(1-xi);
N2 = @(xi) 1/2*(1+xi);
N = @(xi)[N1(xi) N2(xi)];
B=(1/Le)*[-1. 1.];
% gauss quadrature: 2 gauss points
xi = [-1/sqrt(3) 1/sqrt(3)];
w = [1 1];
J=Le/2;

M = rho*A*( (w(1)*N(xi(1))'*N(xi(1))*J)+ (w(2)*N(xi(2))'*N(xi(2))*J) ); % mass matrix
K = E*(w(1)*B'*B*J + w(2)*B'*B*J);

% K(1,1) = 0; K(1,2) = 0; K(2,1) = 0; % encastred beam  
eigval = eig(inv(M)*K);
omega= sqrt(eigval(2));
h=[0.001:0.001:10];
Omega = sqrt(eigval(2))*h;
l_Omega = Omega - 4.0 ;
[val_min,i_min] = min(abs(l_Omega)) ;

for i=1:length(h)
    A = [eye(size(M^-1*K)) gamma*h(i)*M^-1*K;
         zeros(size(M^-1*K)) eye(size(M^-1*K))+h(i)^2*beta*M^-1*K];
     B = -[-eye(size(M^-1*K)) (1-gamma)*h(i)*M^-1*K;
         -eye(size(M^-1*K))*h(i) -eye(size(M^-1*K))+h(i)^2*(1/2-beta)*M^-1*K];
    % Calculation of the amplification matrix
    A_xi(:,:,i) = inv(A) * B;
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
        plot(eigen_val((i-1)*2+1,3)...
            ,eigen_val((i-1)*2+2,3),'r+');
         plot(eigen_val((i-1)*2+1,4)...
             ,eigen_val((i-1)*2+2,4),'og');
        hold on ;
    end
end
title('Eigen values in imaginary plan')

grid on ;
axis equal ;
axis([-2 2 -2 2]) ;

% spectral radius
figure(2) ;
clf ;
for i = 1:length(h)
    if mod(i,20)==0 
        semilogx(Omega(i),lambda_xi(i,3),'+r') ;
        semilogx(Omega(i),lambda_xi(i,4),'og') ;
        hold on ;
    end
end
title('Spectral radius')

grid on ;

%% Periodicity error
figure(3) ;
clf ;
for i = 1:i_min
    plot(Omega(i),PE_xi(i,3),'+r') ;
    plot(Omega(i),PE_xi(i,4),'og') ;
    hold on ;
end
title('Relative periodicity error')

grid on ;


% numerical damping ratio
figure(4) ;
clf ;
for i = 1:i_min
    plot(Omega(i),AD_xi(i,3),'+r') ;
    plot(Omega(i),AD_xi(i,4),'og') ;
    hold on ;
end
title('Numerical Damping ratio')

grid on ;

