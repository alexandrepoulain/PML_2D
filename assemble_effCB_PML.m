function [CCB] = assemble_effCB_PML(bulk,mu,cs,b_scal,dt,XB,TB,a0,b0,x0,L2,h,n,dof)
% This function assembles the effective damping matrix for the 2D PML
% Inputs: k: bulk modulus
%         mu: shear coefficient (Lam�)
%         cs: shear wave velocity
%         b_scal: scalar
%         dt: time step
%         XB: coordinates of the nodes in PML
%         TB: connectivity in PML
%         a0: coefficient of attenuation for evanescent waves
%         b0: coefficient of attenuation for propagating waves
%         x0: beginning of the PML
%         L2: length of the PML in x direction
%         h: length of the medium in y direction
%         n: order of the attenuation functions
%         dof: deggres of freedom for a node
% plane strain
D  = [bulk+4/3*mu bulk-2*mu/3 0
	  bulk-2*mu/3 bulk+4/3*mu 0
	    0      0        mu]; 
D = 1/dt*D;
% 3 gauss  points in each direction
ng = 3;
w = [5/9*5/9 5/9*8/9 5/9*5/9 8/9*5/9 8/9*8/9 8/9*5/9 5/9*5/9 5/9*8/9 5/9*5/9];
eta = [-sqrt(3/5) -sqrt(3/5) -sqrt(3/5) 0.         0. 0.        sqrt(3/5)  sqrt(3/5) sqrt(3/5)];
ksi = [-sqrt(3/5) 0.         sqrt(3/5)  -sqrt(3/5) 0. sqrt(3/5) -sqrt(3/5) 0.        sqrt(3/5)];
% material derivative
dN = @(eta,ksi)1/4*[-(1-eta)  (1-eta)  (1+eta) -(1+eta); 
	                -(1-ksi) -(1+ksi)  (1+ksi)  (1-ksi)];
% rotation matrix
theta = 0;
Q = [cos(theta) sin(theta); -sin(theta) cos(theta)];
% Function of attenuation in 2D
fe1 = @(x) a0(1)*((x-x0)/L2)^n * (x>x0); 
fe2 = @(y) a0(2)*((y-h)/L2)^n * (y>h);
fp1 = @(x) b0(1)*((x-x0)/L2)^n*(x>x0); 
fp2 = @(y) b0(2)*((y-h)/L2)^n*(y>h);
% stability
% fe1 = @(x) a0(1); 
% fe2 = @(y) a0(2);
% fp1 = @(x) b0(1); 
% fp2 = @(y) b0(2);
% Matrix Fp, Fe, Fte, Ftp
Fe = @(x) Q'*[1+fe1(x(1)) 0; 0 1+fe2(x(2))]*Q;
Fp = @(x) Q'*[(cs/b_scal)*fp1(x(1)) 0; 0 (cs/b_scal)*fp2(x(2))]*Q;
Fte = @(x) Q'*[1+fe2(x(2)) 0; 0 1+fe1(x(1))]*Q;
Ftp = @(x) Q'*[(cs/b_scal)*fp2(x(2)) 0; 0 (cs/b_scal)*fp1(x(1))]*Q;
% Matrix Fl, Feps
Fl = @(x)(Fp(x)+(Fe(x)/dt))^(-1);
Feps = @(x) Fe(x)*Fl(x);
% number of nodes in an element
nn = size(TB,2)-1;
% For each element in the PML (each row in TB)
CCB = zeros(size(XB,1)*dof,size(XB,1)*dof);
for ii = 1:size(TB,1)
   % Get back the position of the nodes
   Xe = XB(TB(ii,1:nn),:);
   cce = zeros(dof*nn,dof*nn);
   % For each node
   Bt = zeros(3,8);Beps = zeros(3,8);
   for k=1:ng^2
       M = (1/4)*[((1-ksi(k))*(1-eta(k))) ((1+ksi(k))*(1-eta(k))) ((1+ksi(k))*(1+eta(k))) ((1-ksi(k))*(1+eta(k)))] ;
       pos = M*Xe; % position of the node in natural coordinates 
       % calculation of the matrices
       tempFl = Fl(pos);
       tempdN = dN(eta(k),ksi(k));
       J = tempdN*Xe;
       
       tempdN = J\tempdN;
       tempFesp = Feps(pos);
       tempFte = Fte(pos); tempFtp = Ftp(pos);
       % calculation of Nl,Bq,Nte and Ntp
       ind = 1;
       for kk=1:4
           Nl = [tempFl(1,1)*tempdN(1,kk)+tempFl(1,2)*tempdN(2,kk);
                 tempFl(2,1)*tempdN(1,kk)+tempFl(2,2)*tempdN(2,kk)];
           Beps(:,ind:ind+1) = [tempFesp(1,1)*Nl(1) tempFesp(2,1)*Nl(1);
                 tempFesp(1,2)*Nl(2) tempFesp(2,2)*Nl(2);
                 tempFesp(1,1)*Nl(2)+tempFesp(1,2)*Nl(1) tempFesp(2,1)*Nl(2)+tempFesp(2,2)*Nl(1)];  
           Nte = [tempFte(1,1)*tempdN(1,kk)+tempFte(2,1)*tempdN(2,kk);
                  tempFte(1,2)*tempdN(1,kk)+tempFte(2,2)*tempdN(2,kk)];
           Ntp = [tempFtp(1,1)*tempdN(1,kk)+tempFtp(2,1)*tempdN(2,kk);
                  tempFtp(1,2)*tempdN(1,kk)+tempFtp(2,2)*tempdN(2,kk)]; 
           % Calculation of B tilde matrices
           Bte = [Nte(1) 0; 0 Nte(2); Nte(2) Nte(1)];
           Btp = [Ntp(1) 0; 0 Ntp(2); Ntp(2) Ntp(1)];
           Bt(:,ind:ind+1) = Bte + dt*Btp;
           ind=ind+2;
       end
       % for each gauss point
       J = dN(eta(k),ksi(k))*Xe;
       cce = cce + w(k)*(Bt'*D*Beps)*abs(det(J));       
   end

   % assemble
   % Define global address vector ig( ) for element dofs.
   Te = TB(ii,:);
   ig  = zeros(1,nn*dof);
   for i = 1:nn
     for j = 1:dof
       ig((i-1)*dof+j) = (Te(i)-1)*dof + j;
     end
   end
   % Add element contributions.
   for i = 1:nn*dof
    for j = 1:nn*dof
       CCB(ig(i),ig(j)) = CCB(ig(i),ig(j)) + cce(i,j);
     end
   end
end
end

