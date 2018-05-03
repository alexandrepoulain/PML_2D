function [Mat_B, Mat_Bp] = fun_B(a0,b0,L2,n,x0,h,cs,b_scal,Xe,eta,ksi,dt,k)
% This function calculates the matrices B and Bp
% Inputs: a0,b0 : coefficients of attenuation in both directions
%         x0: start of PML in x direction
%         y0: start of PML in y direction
%         L2: Length of PML (in x direction)
%         h:  Length of PML (in y direction)
%         n: order of attenuation on PML
%         cs: celerity of S-wave
%         b_scal: scalar
%         Xe: position of the nodes in the concerned element
%         eta, ksi: gauss points position
%         k: at which gauss point are we in the element
%         dt: time step on PML

% angle of incidence
theta = 0;
% Rotation matrix
Q = [cos(theta) sin(theta);
    -sin(theta) cos(theta)];
% material derivative
dN = @(k)1/4*[-(1-eta(k))  (1-eta(k))  (1+eta(k)) -(1+eta(k)); 
	          -(1-ksi(k)) -(1+ksi(k))  (1+ksi(k))  (1-ksi(k))];
% Function of attenuation in 2D
% fe1 = @(x) a0(1)*((x-x0)/L2)^n * (x>x0); 
% fe2 = @(y) a0(2)*((y-h)/L2)^n * (y>h);
% fp1 = @(x) b0(1)*((x-x0)/L2)^n*(x>x0); 
% fp2 = @(y) b0(2)*((y-h)/L2)^n*(y>h);
% stability
fe1 = @(x) a0(1); 
fe2 = @(y) a0(2);
fp1 = @(x) b0(1); 
fp2 = @(y) b0(2);
% Matrix Fp, Fe, Fte, Ftp
Fe = @(x) Q'*[1+fe1(x(1)) 0; 0 1+fe2(x(2))]*Q;
Fp = @(x) Q'*[(cs/b_scal)*fp1(x(1)) 0; 0 (cs/b_scal)*fp2(x(2))]*Q;
Fte = @(x) Q'*[1+fe2(x(2)) 0; 0 1+fe1(x(1))]*Q;
Ftp = @(x) Q'*[(cs/b_scal)*fp2(x(2)) 0; 0 (cs/b_scal)*fp1(x(1))]*Q;
% Matrix Fl, Fq
Fl = @(x)(Fp(x)+(Fe(x)/dt))^(-1);
Feps = @(x) Fe(x)*Fl(x);
Fq = @(x) Fp(x)*Fl(x);
% For each node
Mat_B = zeros(3,8);Bq = zeros(3,8);Mat_Bp = zeros(3,8); 

M = (1/4)*[((1-ksi(k))*(1-eta(k))) ((1+ksi(k))*(1-eta(k))) ((1+ksi(k))*(1+eta(k))) ((1-ksi(k))*(1+eta(k)))] ;
pos = M*Xe; % position of the node in natural coordinates 
% calculation of the matrices
tempFl = Fl(pos);
tempdN = dN(k);
J = tempdN*Xe;
tempdN = J\tempdN;
tempFq = Fq(pos);
tempFte = Fte(pos); tempFtp = Ftp(pos);
ind = 1;
for kk=1:4
   % calculation of Nl,Bq,Nte and Ntp
   Nl = [tempFl(1,1)*tempdN(1,kk)+tempFl(1,2)*tempdN(2,kk);
         tempFl(2,1)*tempdN(1,kk)+tempFl(2,2)*tempdN(2,kk)];
   Bq(:,ind:ind+1) = [tempFq(1,1)*Nl(1) tempFq(2,1)*Nl(1);
         tempFq(1,2)*Nl(2) tempFq(2,2)*Nl(2);
         tempFq(1,1)*Nl(2)+tempFq(1,2)*Nl(1) tempFq(2,1)*Nl(2)+tempFq(2,2)*Nl(1)];    
   Nte = [tempFte(1,1)*tempdN(1,kk)+tempFte(2,1)*tempdN(2,kk);
          tempFte(1,2)*tempdN(1,kk)+tempFte(2,2)*tempdN(2,kk)];
   Ntp = [tempFtp(1,1)*tempdN(1,kk)+tempFtp(2,1)*tempdN(2,kk);
          tempFtp(1,2)*tempdN(1,kk)+tempFtp(2,2)*tempdN(2,kk)];   
   % Calculation of B tilde matrices
   Bte = [Nte(1) 0; 0 Nte(2); Nte(2) Nte(1)];
   Btp = [Ntp(1) 0; 0 Ntp(2); Ntp(2) Ntp(1)];
   Mat_B(:,ind:ind+1) = Bte + dt*Btp;
   Mat_Bp(:,ind:ind+1) = Btp;

   ind=ind+2;
end

end

