function [P] = internal_forces_PML(XB,TB,sig,Sig,dof)
% This function calculate the vector of internal forces in the PML.
% Inputs: sig: Matrix of the stress at each Gauss point 
%         Sig: Matrix of the integral of the stress at each gauss point
%         XB: Coordinates of each nodes in the PML
%         TB: Connectivity matrix
%         dof: Number of degree of freedom


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
Feps = @(x) Fe(x)*Fl(x);
% number of nodes in an element
nn = size(TB,2)-1;
% For each element in the PML (each row in TB)
P = zeros(size(XB,1)*dof,1);
for ii = 1:size(TB,1)
   % Get back the position of the nodes
   Xe = XB(TB(ii,1:nn),:);
   sig_e(:,1) = sig(ii*dof:ii*dof+dof*nn,1); % stress at each gauss point
   Sig_e(:,1) = Sig(ii*dof:ii*dof+dof*nn,1); % stress at each gauss point
   pe = zeros(dof*nn,1);
   % For each node
   Bt = zeros(3,8);Beps = zeros(3,8);
   ind=1;
   for k=1:nn
       M = (1/4)*[((1-ksi(k))*(1-eta(k))) ((1+ksi(k))*(1-eta(k))) ((1+ksi(k))*(1+eta(k))) ((1-ksi(k))*(1+eta(k)))] ;
       pos = M*Xe; % position of the node in natural coordinates 
       % calculation of the matrices
       tempFl = Fl(pos);
       tempdN = dN(eta,ksi,k);
       J = tempdN*Xe;
       tempdN = J\tempdN;
       tempFesp = Feps(pos);
       tempFte = Fte(pos); tempFtp = Ftp(pos);
       % calculation of Nl,Bq,Nte and Ntp
       Nl = [tempFl(1,1)*tempdN(1,k)+tempFl(1,2)*tempdN(2,k);
             tempFl(2,1)*tempdN(1,k)+tempFl(2,2)*tempdN(2,k)];
       Beps(:,ind:ind+1) = [tempFesp(1,1)*Nl(1) tempFesp(2,1)*Nl(1);
             tempFesp(1,2)*Nl(2) tempFesp(2,2)*Nl(2);
             tempFesp(1,1)*Nl(2)+tempFesp(1,2)*Nl(1) tempFesp(2,1)*Nl(2)+tempFesp(2,2)*Nl(1)];    
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
       pe = pe + w(k)*(Bt'*sig_e + Btp'*Sig_e)*abs(det(J));
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
   for i = 1:nn
    for j = 1:dof
       P(ig(i)+ig(j),1) = P(ig(i)+ig(j),1) + pe((i-1)*dof+j);
     end
   end

end

end

