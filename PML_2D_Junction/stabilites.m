% domaine de stabilit� des sch�mas temporels 
clear all;
close all;
%% Space parameters
Le = 1;
rho = 1700; % mass density
mu_e = rho;
S = 1;
Lpml = 50; % length of the PML
E=10e5; % Young modulus
cs = Lpml;
xi1 = -1/sqrt(3);
xi2 = 1/sqrt(3);
N1 = @(xi) 1/2*(1-xi);
N2 = @(xi) 1/2*(1+xi);
w1 =1;
w2 =1;
R=(1/Le)*[-1. 1];
%% Assembling matrices for medium
% assembling M matrix
M = zeros(2,2);
M = mu_e*Le/2*(w1*[N1(xi1)^2 N1(xi1)*N2(xi1) ; N1(xi1)*N2(xi1) N2(xi1)^2] + ...
    w2*[N1(xi2)^2 N1(xi2)*N2(xi2); N1(xi2)*N2(xi2) N2(xi2)^2]);
% assembling K matrix
K = zeros(2,2);
% submatrix
K = E * Le/2 *(w1*(R)'*R +w2*(R)'*R);
H = inv(M)*K;
eigenv = eig(H);

% matrice d'amplification dans le cas Newmark non amorti
% param�tres de Newmark
gamma=0.5;beta=0.25;
omega = eigenv(2);
h=[0.001:0.001:10];
dt = h;
Omega = omega*h ;
%l_Omega = Omega - 2.5 ;
l_Omega = Omega - 4.0 ;
[val_min,i_min] = min(abs(l_Omega)) ;

%% assembling little matrices
% Shape functions in normed space


fp = 10;

M = rho*Le/2*S*([N1(xi1)^2 N1(xi1)*N2(xi1) ; N1(xi1)*N2(xi1) N2(xi1)^2] ...
        + [N1(xi2)^2 N1(xi2)*N2(xi2); N1(xi2)*N2(xi2) N2(xi2)^2]);

C = rho/Lpml *cs *S* Le/2 * (fp.*[N1(xi1)^2 N1(xi1)*N2(xi1) ; N1(xi1)*N2(xi1) N2(xi1)^2] ...
        + fp.*[N1(xi2)^2 N1(xi2)*N2(xi2); N1(xi2)*N2(xi2) N2(xi2)^2]);



%% 
for i=1:length(h)
    alpha = 1+cs/Lpml*h(i)*fp ;

    K = E* Le/2/alpha .*((R)'*R +(R)'*R);
    %%%%% ------------- cas amorti Newmark dans le cas gamma=1/2 et beta=1/4
    tempB = 1/Le.*[-1 1; -1 1];
    Ap = 1/2*cs/Lpml*E*S*fp/alpha.*[-1 -1; 1 1];
    A = [ eye(2)+gamma*h(i).*M^-1*C gamma*h(i).*M^-1*K -gamma*h(i).*M^-1*Ap;
        beta*h(i)^2.*M^-1*C eye(2)+beta*h(i)^2.*M^-1*K -beta*h(i)^2.*M^-1*Ap;
        zeros(2,2) -h(i)/alpha.*tempB eye(2)];

    B = -[ -eye(2)+(1-gamma)*h(i).*M^-1*C (1-gamma)*h(i).*M^-1*K  -(1-gamma)*h(i).*M^-1*Ap;
        -h(i).*eye(2)+(1/2-beta)*h(i)^2.*M^-1*C -eye(2)+(1/2-beta)*h(i)^2.*M^-1*K  -(1/2-beta)*h(i)^2.*M^-1*Ap;
      zeros(2,2) zeros(2,2) -eye(2)+h(i)*cs/Lpml*fp/alpha.*eye(2)];
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
        %plot(eigen_val((i-1)*2+1,1)...
        %    ,eigen_val((i-1)*2+2,1),'+r');
        %plot(eigen_val((i-1)*2+1,2)...
        %    ,eigen_val((i-1)*2+2,2),'og');
        hold on ;
    end
end
grid on ;
axis equal ;
axis([-2 2 -2 2]) ;

% spectral radius
figure(2) ;
clf ;
for i = 1:length(dt)
    if mod(i,20)==0 
        semilogx(Omega(i),lambda_xi(i,:),'+') ;
        %semilogx(Omega(i),lambda_xi(i,1),'+r') ;
        %semilogx(Omega(i),lambda_xi(i,2),'og') ;

        hold on ;
    end
end
grid on ;

%% Periodicity error
figure(3) ;
clf ;
for i = 1:i_min
    plot(Omega(i),PE_xi(i,:),'+') ;
    %plot(Omega(i),PE_xi(i,1),'+r') ;
    %plot(Omega(i),PE_xi(i,2),'og') ;
    hold on ;
end
grid on ;


% numerical damping ratio
figure(4) ;
clf ;
for i = 1:i_min
    plot(Omega(i),AD_xi(i,:),'+') ;
    %plot(Omega(i),AD_xi(i,1),'+r') ;
    %plot(Omega(i),AD_xi(i,2),'og') ;
    hold on ;
end
grid on ;
