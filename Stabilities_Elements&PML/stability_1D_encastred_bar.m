%% 1 dimentional SDOF bar element Stability

clear all; close all;
% Newmark parameters
gamma = 1/2; beta = 1/4; % Implicit
gamma = 1/2; beta = 0;   % Explicit
rho = 1700; % mass density
E = 10e3;
omega=E/rho;
h=[0.001:0.001:10];
Omega = omega*h;
l_Omega = Omega - 4.0 ;
[val_min,i_min] = min(abs(l_Omega)) ;

for i=1:length(h)
    A = [1 omega^2*gamma*h(i);
         0 h(i)^2*beta*omega^2+1];
     B = -[-1 omega^2*(1-gamma)*h(i);
         -h(i) -1+h(i)^2*(1/2-beta)*omega^2];
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
       PE_xi(i,k) = ( omega * h(i) / abs(angle(D(k,i))) ) - 1;
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
        plot(eigen_val((i-1)*2+1,1)...
            ,eigen_val((i-1)*2+2,1),'+r');
        plot(eigen_val((i-1)*2+1,2)...
            ,eigen_val((i-1)*2+2,2),'og');
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
        semilogx(Omega(i),lambda_xi(i,1),'+r') ;
        semilogx(Omega(i),lambda_xi(i,2),'og') ;
        hold on ;
    end
end
title('Spectral radius')

grid on ;

%% Periodicity error
figure(3) ;
clf ;
for i = 1:i_min
    plot(Omega(i),PE_xi(i,1),'+r') ;
    plot(Omega(i),PE_xi(i,2),'og') ;
    hold on ;
end
title('Relative periodicity error')

grid on ;


% numerical damping ratio
figure(4) ;
clf ;
for i = 1:i_min
    plot(Omega(i),AD_xi(i,1),'+r') ;
    plot(Omega(i),AD_xi(i,2),'og') ;
    hold on ;
end
title('Numerical Damping ratio')

grid on ;

