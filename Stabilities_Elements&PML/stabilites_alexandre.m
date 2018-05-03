% domaine de stabilite des sch�mas temporels 
clear all;
close all;
% matrice d'amplification dans le cas Newmark non amorti
% parametres de Newmark
%gamma=0.5;beta=0.25;
% param�tres de Newmark
alpha_Newmark = 0.0 ;
gamma = 1/2 + alpha_Newmark ;
beta = 1/4*(1+alpha_Newmark)^2 ;
% cas explicite 
gamma = 1/2 ; beta = 0 ;

%omega=2*pi;
%h=[0.001:0.001:10];
% h=[[0.00001:0.00001:0.0001] [0.002:0.001:10]] ;
% Omega = omega*h ;
% %l_Omega = Omega - 2.5 ;
% l_Omega = Omega - 4.0 ;
% [val_min,i_min] = min(abs(l_Omega)) ;

%% assembling little matrices
% Shape functions in normed space
Le = 1;
rho = 1700; % mass density
S = 1;
Lpml = 50; % length of the PML
E=10e5; % Young modulus
cs = Lpml;
xi1 = -1/sqrt(3);
xi2 = 1/sqrt(3);
N1 = @(xi) 1/2*(1-xi);
N2 = @(xi) 1/2*(1+xi);
R=(1/Le)*[-1. 1];
fp = 15;

M = rho*Le/2*S*([N1(xi1)^2 N1(xi1)*N2(xi1) ; N1(xi1)*N2(xi1) N2(xi1)^2] ...
        + [N1(xi2)^2 N1(xi2)*N2(xi2); N1(xi2)*N2(xi2) N2(xi2)^2]);

C = rho/Lpml *cs *S* Le/2 * (fp.*[N1(xi1)^2 N1(xi1)*N2(xi1) ; N1(xi1)*N2(xi1) N2(xi1)^2] ...
        + fp.*[N1(xi2)^2 N1(xi2)*N2(xi2); N1(xi2)*N2(xi2) N2(xi2)^2]);


% cas Newmark sans amortissement 
% matrice de rigidit� 
B = (1/Le)*[-1. 1] ; 
K0 = E* Le/2 .*((B)'*B +(B)'*B);

[V0,D0] = eig(K0) ;
% V0(:,2) est la signature des modes propres conjugu�s
% taille de la matrice A(h) 
dim = 2 * length(B) ;
K0_rebuild = zeros(2,2);
for k=1:2 
    K0_rebuild(:,:) = K0_rebuild(:,:) +  D0(k,k)*V0(:,k)*V0(:,k)';
end ;

% fr�quence propre (1 mode de corps rigide avec une fr�quence � 0 puis
% l'autre mode)

[V1,D1] = eig(inv(M)*K0) ;
omega=sqrt(D1(2,2)); % 84.0168 rad/s
% essai reconstruction 
MK_rebuild = zeros(2,2);
for k=1:2 
    MK_rebuild(:,:) = MK_rebuild(:,:) +  D1(k,k)*V1(:,k)*V1(:,k)';
end ;
% verifier par rappirt � inv(M)*K0
%h=[0.001:0.001:10];
h=[[0.0001:0.0001:0.1]] ;
Omega = omega*h ;
l_Omega = Omega - 4.0 ;
[val_min,i_min] = min(abs(l_Omega)) ;

C0 = zeros(2,2) ;

% % matrice d'amorissement Rayleigh
% xi = 0.05 ;
% % uniquement fonction de la matrice de Masse
% %alpha_M = 2 * xi * omega ;
% %C0 = alpha_M * M ;
% % % uniquement fonction de la matrice de raideur
% % %alpha_K = 2 * xi / omega ;
% % %C0 = alpha_K * K0 ;
% % % fonction de la matrice de raideur et de masse
% alpha_M = xi * omega ;
% alpha_K = xi / omega ;
% C0 = alpha_M * M + alpha_K * K0 ;

epsilon = 1e-2 ;
for i=1:length(h)
  %h(1) = 0 ;  
  % matrices H1 et H0 : vitesses aux 2 noeuds puis d�placements aux deux noeuds 
  H1 = [ M+gamma*h(i)*C0    gamma*h(i)*K0 ; h(i)^2*beta*C0 (M+h(i)^2*beta*K0) ] ;
  H0 = [ M-h(i)*(1-gamma)*C0    -h(i)*(1-gamma)*K0 ; (h(i)*M - h(i)^2*(1/2-beta)*C0) (M-h(i)^2*(1/2-beta)*K0) ] ;
  % matrice A(h)
  A = inv(H1) * H0 ;
  % test du eigenshuffle
  %[V,D] = eig(A) ;
  %d_newmarkI(i,:) = diag(D) ;
  %
  [V,D] = eigenshuffle(A) ;
  d_newmarkI(i,:) = D ;
  %
  V_newmarkI(:,:,i) = V ;
  % reconnaissance des modes propres relatifs � la d�formation de traction
  % et compression
%   compt = 0 ;
%   for k=1:dim
%       if (abs(V(1,k)+V0(:,2))<epsilon) 
%           VP_P(i,compt+1) = d_newmarkI(i,k);
%           compt = compt + 1 ;
%       end;
%   end;
  %[d_newmarkI(i,:),vect_newmarkI(i,:)]_=eig(A) ;
  
  % abscisse = partie r�elle des valuers propres de la matrice
  % d'amplification A 
  SR(i) = max(abs(d_newmarkI(i,:))) ;
  A_rebuild=zeros(4,4) ;
  % boucle sur les VPs
  for k=1:dim
    % partie r�elle et imaginaire de la VP
    x(i,k) = real(d_newmarkI(i,k)) ;
    y(i,k) = imag(d_newmarkI(i,k)) ;
    % modules spectraux 
    lambda(i,k) =  abs(d_newmarkI(i,k)) ;
    % erreur de p�riodicit�
    PE(i,k) = ( omega * h(i) / abs(angle(d_newmarkI(i,k)) ) ) - 1 ;
    % amortissement num�rique
    AD(i,k)= -log(lambda(i,k)) / abs(angle(d_newmarkI(i,k))) ;
    % verification de la decomposition spectrale 
    A_rebuild = A_rebuild + d_newmarkI(i,k) * V_newmarkI(:,k,i) * V_newmarkI(:,k,i)' ;
  end;
  % plus simplement : A_rebuild = V * diag(D) * V'
end;


% trac� dans un plan partie imaginaire et partie r�elle 
figure(1) ;
clf ;
% for k=1:dim
% plot(x(1:length(x),k),y(1:length(y),k),'+') ;
% hold on ;
% end;
plot(x(1:length(x),1),y(1:length(y),1),'r+') ;
hold on ;
plot(x(1:length(x),2),y(1:length(y),2),'bo') ;
hold on ;
plot(x(1:length(x),3),y(1:length(y),3),'kx') ;
hold on ;
plot(x(1:length(x),4),y(1:length(y),4),'g^') ;
grid on ;
title('Eigenvalues in imaginary plan')
axis equal ;
axis([-2 2 -2 2]) ;

% trac� du rayon spectral 
figure(2) ;
clf ;
semilogx(Omega(1:1:length(Omega)),lambda(1:1:length(Omega),1),'r+') ;
hold on ;
semilogx(Omega(1:1:length(Omega)),lambda(1:1:length(Omega),2),'bo') ;
hold on ;
semilogx(Omega(1:1:length(Omega)),lambda(1:1:length(Omega),3),'kx') ;
hold on ;
semilogx(Omega(1:1:length(Omega)),lambda(1:1:length(Omega),4),'g^') ;
grid on ;
title('Spectral radius')

axis([0.01 max(Omega) 0.85 1.2]) ;


% trac� de l'erreur de p�riodicit�
figure(3) ;
clf ;
%plot(Omega(1:1:i_min),PE(1:1:i_min,1),'r+') ;
%hold on ;
%plot(Omega(1:1:i_min),PE(1:1:i_min,2),'bo') ;
%hold on ;
plot(Omega(1:1:i_min),PE(1:1:i_min,3),'kx') ;
hold on ;
plot(Omega(1:1:i_min),PE(1:1:i_min,4),'g^') ;
title('Relative periodicity error')

grid on ;


% trac� de l'amortissement num�rique
figure(4) ;
clf ;
%plot(Omega(1:1:i_min),AD(1:1:i_min,1),'r+') ;
%hold on ;
%plot(Omega(1:1:i_min),AD(1:1:i_min,2),'bo') ;
%hold on ;
plot(Omega(1:1:i_min),AD(1:1:i_min,3),'kx') ;
hold on ;
plot(Omega(1:1:i_min),AD(1:1:i_min,4),'g^') ;
hold on ;
title('Numerical damping ratio')

grid on ;
disp('sch�ma CAA Newmark')
pente_xi = (log(AD(2,4))-log(AD(1,4))) / (log(Omega(2))-log(Omega(1)))
pente_T = (log(PE(2,4))-log(PE(1,4))) / (log(Omega(2))-log(Omega(1)))



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cas de la PML 1D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

  B = [ -eye(2)+(1-gamma)*h(i).*M^-1*C (1-gamma)*h(i).*M^-1*K  -(1-gamma)*h(i).*M^-1*Ap;
        -h(i).*eye(2)+(1/2-beta)*h(i)^2.*M^-1*C -eye(2)+(1/2-beta)*h(i)^2.*M^-1*K  -(1/2-beta)*h(i)^2.*M^-1*Ap;
      zeros(2,2) zeros(2,2) -eye(2)+h(i)*cs/Lpml*fp/alpha.*eye(2)];

  A_PML = -inv(A) * B ;
  %d_xi(i,:)=eig(A_xi) ;
  [V,D] = eigenshuffle(A_PML) ;
  d_newmarkI_PML(i,:) = D ;
  V_newmarkI_PML(:,:,i) = V ;
  % rayon spectral de la matrice
  SR_PML(i) = max(abs(d_newmarkI_PML(i,:))) ;
  A_rebuild = zeros(6,6) ;
  % boucle sur les VPs
  for k=1:6
    % partie r�elle et imaginaire de la VP
    x_PML(i,k) = real(d_newmarkI_PML(i,k)) ;
    y_PML(i,k) = imag(d_newmarkI_PML(i,k)) ;
    % modules spectraux 
    lambda_PML(i,k) =  abs(d_newmarkI_PML(i,k)) ;
    % erreur de p�riodicit�
    PE_PML(i,k) = ( omega * h(i) / abs(angle(d_newmarkI_PML(i,k)) ) ) - 1 ;
    % amortissement num�rique
    AD_PML(i,k)= -log(lambda_PML(i,k)) / abs(angle(d_newmarkI_PML(i,k))) ;
    % verification de la decomposition spectrale 
    A_rebuild(:,:) = A_rebuild(:,:) + d_newmarkI_PML(i,k) * V_newmarkI_PML(:,k,i) * V_newmarkI_PML(:,k,i)' ;
  end;
end;


% trac� dans un plan partie imaginaire et partie r�elle 
figure(11) ;
clf ;
plot(x_PML(1:length(x_PML),1),y_PML(1:length(y_PML),1),'r+') ;
hold on ;
plot(x_PML(1:length(x_PML),2),y_PML(1:length(y_PML),2),'bo') ;
hold on ;
plot(x_PML(1:length(x_PML),3),y_PML(1:length(y_PML),3),'kx') ;
hold on ;
plot(x_PML(1:length(x_PML),4),y_PML(1:length(y_PML),4),'g^') ;
hold on ;
plot(x_PML(1:length(x_PML),5),y_PML(1:length(y_PML),5),'mp') ;
hold on ;
plot(x_PML(1:length(x_PML),6),y_PML(1:length(y_PML),6),'yh') ;
grid on ;
title('Eigenvalues in imaginary plan')

axis equal ;
axis([-2 2 -2 2]) ;

% trac� du rayon spectral 
figure(12) ;
clf ;
semilogx(Omega(1:1:length(Omega)),lambda_PML(1:1:length(Omega),1),'r+') ;
hold on ;
semilogx(Omega(1:1:length(Omega)),lambda_PML(1:1:length(Omega),2),'bo') ;
hold on ;
semilogx(Omega(1:1:length(Omega)),lambda_PML(1:1:length(Omega),3),'kx') ;
hold on ;
semilogx(Omega(1:1:length(Omega)),lambda_PML(1:1:length(Omega),4),'g^') ;
hold on ;
semilogx(Omega(1:1:length(Omega)),lambda_PML(1:1:length(Omega),5),'mp') ;
hold on ;
semilogx(Omega(1:1:length(Omega)),lambda_PML(1:1:length(Omega),6),'yh') ;
grid on ;
title('Spectral radius')

axis([0.01 max(Omega) 0.0 1.2]) ;



% trac� de l'erreur de p�riodicit�
figure(13) ;
clf ;
% plot(Omega(1:1:i_min),PE_PML(1:1:i_min,1),'r+') ;
% hold on ;
% plot(Omega(1:1:i_min),PE_PML(1:1:i_min,2),'bo') ;
% hold on ;
% plot(Omega(1:1:i_min),PE_PML(1:1:i_min,3),'kx') ;
% hold on ;
% plot(Omega(1:1:i_min),PE_PML(1:1:i_min,4),'g^') ;
% hold on ;
plot(Omega(1:1:i_min),PE_PML(1:1:i_min,5),'mp') ;
hold on ;
plot(Omega(1:1:i_min),PE_PML(1:1:i_min,6),'yh') ;
title('Relative periodicity error')

grid on ;


% trac� de l'amortissement num�rique
figure(14) ;
clf ;
% plot(Omega(1:1:i_min),AD_PML(1:1:i_min,1),'r+') ;
% hold on ;
% plot(Omega(1:1:i_min),AD_PML(1:1:i_min,2),'bo') ;
% hold on ;
% plot(Omega(1:1:i_min),AD_PML(1:1:i_min,3),'kx') ;
% hold on ;
% plot(Omega(1:1:i_min),AD_PML(1:1:i_min,4),'g^') ;
% hold on ;
plot(Omega(1:1:i_min),AD_PML(1:1:i_min,5),'mp') ;
hold on ;
plot(Omega(1:1:i_min),AD_PML(1:1:i_min,6),'yh') ;
grid on ;
title('Numerical damping ratio')


grid on ;
disp('sch�ma CAA Newmark PML')
pente_xi = (log(AD_PML(2,5))-log(AD_PML(1,5))) / (log(Omega(2))-log(Omega(1)))
pente_T = (log(PE_PML(2,5))-log(PE_PML(1,5))) / (log(Omega(2))-log(Omega(1)))

pouet

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r�sultats SDOF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------- Newmark avec amortissement visqueux
figure(1) ;
clf ;
plot(x1_xi(1:20:length(x1_xi)),y1_xi(1:20:length(y1_xi)),'r+') ;
hold on ; 
plot(x2_xi(1:20:length(x2_xi)),y2_xi(1:20:length(y2_xi)),'bo') ;
hold on ; 
plot(x3_xi(1:20:length(x3_xi)),y3_xi(1:20:length(y3_xi)),'go') ;
hold on ; 
plot(x4_xi(1:20:length(x4_xi)),y4_xi(1:20:length(y4_xi)),'ko') ;
hold on ; 
plot(x5_xi(1:20:length(x5_xi)),y5_xi(1:20:length(y5_xi)),'co') ;
hold on ; 
plot(x6_xi(1:20:length(x6_xi)),y6_xi(1:20:length(y6_xi)),'mo') ;
grid on ;
axis equal ;
axis([-2 2 -2 2]) ;

% trac� du rayon spectral 
figure(2) ;
clf ;
semilogx(Omega(1:20:length(Omega)),lambda1_xi(1:20:length(lambda1_xi)),'r+') ;
hold on ;
semilogx(Omega(1:20:length(Omega)),lambda2_xi(1:20:length(lambda2_xi)),'bo') ;
hold on ;
semilogx(Omega(1:20:length(Omega)),lambda3_xi(1:20:length(lambda3_xi)),'go') ;
hold on ;
semilogx(Omega(1:20:length(Omega)),lambda4_xi(1:20:length(lambda4_xi)),'ko') ;
hold on ;
semilogx(Omega(1:20:length(Omega)),lambda5_xi(1:20:length(lambda5_xi)),'co') ;
hold on ;
semilogx(Omega(1:20:length(Omega)),lambda6_xi(1:20:length(lambda6_xi)),'mo') ;
grid on ;


% trac� de l'erreur de p�riodicit�
figure(3) ;
clf ;
plot(Omega(1:20:i_min),PE1_xi(1:20:i_min),'r+') ;
hold on ;
plot(Omega(1:20:i_min),PE2_xi(1:20:i_min),'bo') ;
hold on ;
plot(Omega(1:20:i_min),PE3_xi(1:20:i_min),'go') ;
hold on ;
plot(Omega(1:20:i_min),PE4_xi(1:20:i_min),'ko') ;
hold on ;
plot(Omega(1:20:i_min),PE5_xi(1:20:i_min),'co') ;
hold on ;
plot(Omega(1:20:i_min),PE6_xi(1:20:i_min),'mo') ;
grid on ;


% trac� de l'amortissement num�rique
figure(4) ;
clf ;
plot(Omega(1:20:i_min),AD1_xi(1:20:i_min),'r+') ;
hold on ;
plot(Omega(1:20:i_min),AD2_xi(1:20:i_min),'bo') ;
hold on ;
plot(Omega(1:20:i_min),AD3_xi(1:20:i_min),'go') ;
hold on ;
plot(Omega(1:20:i_min),AD4_xi(1:20:i_min),'ko') ;
hold on ;
plot(Omega(1:20:i_min),AD5_xi(1:20:i_min),'co') ;
hold on ;
plot(Omega(1:20:i_min),AD6_xi(1:20:i_min),'mo') ;

grid on ;
