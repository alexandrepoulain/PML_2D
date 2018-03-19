function [KKB] = assemble_effKB_PML(k,mu,cs,b_scal,dt,XB,TB,a0,b0,x0,L2,h,n,dof)
% This function assembles the effective stiffness matrix for the 2D PML
% Inputs: k: bulk modulus
%         mu: shear coefficient (Lamé)
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
D  = [k+4/3*mu k-2*mu/3 0
	  k-2*mu/3 k+4/3*mu 0
	    0      0        mu]; 

% 3 gauss  points in each direction
ngauss  = 3;
w = [5/9*5/9 5/9*8/9 5/9*5/9 8/9*5/9 8/9*8/9 8/9*5/9 5/9*5/9 5/9*8/9 5/9*5/9];
eta = [-sqrt(3/5) -sqrt(3/5) -sqrt(3/5) 0.         0. 0.        sqrt(3/5)  sqrt(3/5) sqrt(3/5)];
ksi = [-sqrt(3/5) 0.         sqrt(3/5)  -sqrt(3/5) 0. sqrt(3/5) -sqrt(3/5) 0.        sqrt(3/5)];
% material derivative
dN = @(eta,ksi,k)1/4*[-(1-eta(k))  (1-eta(k))  (1+eta(k)) -(1+eta(k)); 
	                -(1-ksi(k)) -(1+ksi(k))  (1+ksi(k))  (1-ksi(k))];
% Function of attenuation in 2D
fe1 = @(x) a0(1)*((x-x0)/L2)^n; 
fe2 = @(y) a0(2)*((y-h)/L2)^n;
fp1 = @(x) b0(1)*((x-x0)/L2)^n; 
fp2 = @(y) b0(2)*((y-h)/L2)^n;
% Matrix Fp, Fe, Fte, Ftp
Fe = @(x) [1+fe1(x(1)) 0; 0 1+fe2(x(2))];
Fp = @(x) [(cs/b_scal)*fp1(x(1)) 0; 0 (cs/b_scal)*fp2(x(2))];
Fte = @(x) [1+fe2(x(2)) 0; 0 1+fe1(x(1))];
Ftp = @(x) [(cs/b_scal)*fp2(x(2)) 0; 0 (cs/b_scal)*fp1(x(1))];
% Matrix Fl, Fq
Fl = @(x)(Fp(x)+(Fe(x)/dt))^(-1);
Fq = @(x) Fp(x)*Fl(x);
% number of nodes in an element
nn = size(TB,2)-1;
% For each element in the PML (each row in TB)
KKB = zeros(size(XB,1)*dof,size(XB,1)*dof);
for ii = 1:size(TB,1)
   % Get back the position of the nodes
   Xe = XB(TB(ii,1:nn),:);
   kke = zeros(dof*nn,dof*nn);
   % For each node
   Bt = zeros(3,8);Bq = zeros(3,8);
   ind=1;
   for k=1:nn
       M = (1/4)*[((1-ksi(k))*(1-eta(k))) ((1+ksi(k))*(1-eta(k))) ((1+ksi(k))*(1+eta(k))) ((1-ksi(k))*(1+eta(k)))] ;
       pos = M*Xe; % position of the node in natural coordinates 
       % calculation of the matrices
       tempFl = Fl(pos);
       tempdN = dN(eta,ksi,k);
       J = tempdN*Xe;
       tempdN = J\tempdN;
       tempFq = Fq(pos);
       tempFte = Fte(pos); tempFtp = Ftp(pos);
       % calculation of Nl,Bq,Nte and Ntp
       Nl = [tempFl(1,1)*tempdN(1,k)+tempFl(1,2)*tempdN(2,k);
             tempFl(2,1)*tempdN(1,k)+tempFl(2,2)*tempdN(2,k)];
       Bq(:,ind:ind+1) = [tempFq(1,1)*Nl(1) tempFq(2,1)*Nl(1);
             tempFq(1,2)*Nl(2) tempFq(2,2)*Nl(2);
             tempFq(1,1)*Nl(2)+tempFq(1,2)*Nl(1) tempFq(2,1)*Nl(2)+tempFq(2,2)*Nl(1)];    
       Nte = [tempFte(1,1)*tempdN(1,k)+tempFte(2,1)*tempdN(2,k);
              tempFte(1,2)*tempdN(1,k)+tempFte(2,2)*tempdN(2,k)];
       Ntp = [tempFtp(1,1)*tempdN(1,k)+tempFtp(2,1)*tempdN(2,k);
              tempFtp(1,2)*tempdN(1,k)+tempFtp(2,2)*tempdN(2,k)];   
       % Calculation of B tilde matrices
       Bte = [Nte(1) 0; 0 Nte(2); Nte(1) Nte(2)];
       Btp = [Ntp(1) 0; 0 Ntp(2); Ntp(1) Ntp(2)];
       Bt(:,ind:ind+1) = Bte + dt*Btp;
       ind=ind+2;
   end
   for k=1:ngauss^2
       % for each gauss point
       J = tempdN*Xe;
       kke = kke + w(k)*1/dt*(Bt'*D*Bq)*abs(det(J));
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
       KKB(ig(i),ig(j)) = KKB(ig(i),ig(j)) + kke(i,j);
     end
   end
end
end

