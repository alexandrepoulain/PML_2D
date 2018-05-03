function [Mass, Stiff, Damp] = construct_global(X,T,G,rho,cs,b_scal,bulk,mu,L2,x0,h,n,dt,dof,ng,a0,b0)
% This function construct the global matrices to bind both sub-domain
% together
% Inputs: X,T,G
% plane strain
nne = size(T,2)-1;
Nef = size(T,1);
D  = [bulk+4/3*mu bulk-2*mu/3 0
	  bulk-2*mu/3 bulk+4/3*mu 0
	    0      0        mu]; 
% Gauss Quadrature 2 points in each direction
switch ng
    case 2
        w = [1 1 1 1];
        eta = [-1/sqrt(3) -1/sqrt(3) 1/sqrt(3) 1/sqrt(3)];
        ksi = [-1/sqrt(3) 1/sqrt(3) 1/sqrt(3) -1/sqrt(3)];
    case 3
        w = [5/9*5/9 5/9*8/9 5/9*5/9 8/9*5/9 8/9*8/9 8/9*5/9 5/9*5/9 5/9*8/9 5/9*5/9];
        eta = [-sqrt(3/5) -sqrt(3/5) -sqrt(3/5) 0.         0. 0.        sqrt(3/5)  sqrt(3/5) sqrt(3/5)];
        ksi = [-sqrt(3/5) 0.         sqrt(3/5)  -sqrt(3/5) 0. sqrt(3/5) -sqrt(3/5) 0.        sqrt(3/5)];
end
% Nodal shape functions
N1 = @(xi, eta) 1/4*(1-xi)*(1-eta);
N2 = @(xi, eta) 1/4*(1+xi)*(1-eta);
N3 = @(xi, eta) 1/4*(1+xi)*(1+eta);
N4 = @(xi, eta) 1/4*(1-xi)*(1+eta);
% Shape functions
N = @(k) [N1(ksi(k),eta(k)) 0 N2(ksi(k),eta(k)) 0 N3(ksi(k),eta(k)) 0 N4(ksi(k), eta(k)) 0;
    0 N1(ksi(k),eta(k)) 0 N2(ksi(k),eta(k)) 0 N3(ksi(k), eta(k)) 0 N4(ksi(k),eta(k))];
% Material derivative
dN = @(k) [ -(1-eta(k))  (1-eta(k))  (1+eta(k)) -(1+eta(k));
	           -(1-ksi(k)) -(1+ksi(k))  (1+ksi(k))  (1-ksi(k)) ]/4;
% Transform matrix 
M = @(k) (1/4)*[((1-ksi(k))*(1-eta(k))) ((1+ksi(k))*(1-eta(k))) ((1+ksi(k))*(1+eta(k))) ((1-ksi(k))*(1+eta(k)))] ;

% rotation matrix
theta = 0;
Q = [cos(theta) sin(theta); -sin(theta) cos(theta)];
% Function of attenuation in 2D
fe1 = @(x) a0(1)*((x-x0)/L2)^n * (x>x0); 
fe2 = @(y) a0(2)*((y-h)/L2)^n * (y>h);
fp1 = @(x) b0(1)*((x-x0)/L2)^n*(x>x0); 
fp2 = @(y) b0(2)*((y-h)/L2)^n*(y>h);
% Matrix Fp, Fe, Fte, Ftp
Fe = @(x) Q'*[1+fe1(x(1)) 0; 0 1+fe2(x(2))]*Q;
Fp = @(x) Q'*[(cs/b_scal)*fp1(x(1)) 0; 0 (cs/b_scal)*fp2(x(2))]*Q;
Fte = @(x) Q'*[1+fe2(x(2)) 0; 0 1+fe1(x(1))]*Q;
Ftp = @(x) Q'*[(cs/b_scal)*fp2(x(2)) 0; 0 (cs/b_scal)*fp1(x(1))]*Q;
% Matrix Fl, Feps
Fl = @(x)(Fp(x)+(Fe(x)/dt))^(-1);
Feps = @(x) Fe(x)*Fl(x);
Fq = @(x) Fp(x)*Fl(x);
fm = @(x) (1+fe1(x(1)))*(1+fe2(x(2)));
fc = @(x) (1+fe1(x(1)))*fp2(x(2))+(1+fe2(x(2)))*fp1(x(1));
fk = @(x) fp1(x(1))*fp2(x(2));

%% Assemble mass matrix 
Mass = zeros(size(X,1)*dof, size(X,1)*dof);
Stiff = zeros(size(X,1)*dof, size(X,1)*dof);
Damp = zeros(size(X,1)*dof, size(X,1)*dof);
% For each element
for i = 1:Nef
   % Position of the nodes
   Xe = X(T(i,1:nne),:);
   Te = T(i,:);
   Ge = G{i};
   % Calculation of submatrix
   me = zeros(2*nne);
   ce = zeros(2*nne);
   ke = zeros(2*nne);
   % test if the element is part of the medium or pml
   switch Ge{1}
       case 'MED'
           for k=1:ng^2
               J = dN(k)*Xe;
               global_dN = J\dN(k);
               B  = [  global_dN(1,1)    0    global_dN(1,2)   0     global_dN(1,3)    0    global_dN(1,4)   0
	               0    global_dN(2,1)    0    global_dN(2,2)    0    global_dN(2,3)    0    global_dN(2,4)
	            global_dN(2,1) global_dN(1,1) global_dN(2,2) global_dN(1,2) global_dN(2,3) global_dN(1,3) global_dN(2,4) global_dN(1,4) ];
               me = me + rho*w(k)*( (N(k))'*N(k) )*abs(det(J));
               ke = ke + w(k)*( B'*D*B )*abs(det(J));
           end

       case 'PML'
           for k=1:ng^2
               pos = M(k)*Xe; % position of the node in natural coordinates 
               % calculation of the matrices
               tempFl = Fl(pos);
               tempdN = dN(k);
               J = tempdN*Xe;
               tempdN = J\tempdN;
               tempFq = Fq(pos); tempFesp = Feps(pos);
               tempFte = Fte(pos); tempFtp = Ftp(pos);
               ind = 1;
               for kk=1:nne
                   % calculation of Nl,Bq,Nte and Ntp
                   Nl = [tempFl(1,1)*tempdN(1,kk)+tempFl(1,2)*tempdN(2,kk);
                         tempFl(2,1)*tempdN(1,kk)+tempFl(2,2)*tempdN(2,kk)];
                   Bq(:,ind:ind+1) = [tempFq(1,1)*Nl(1) tempFq(2,1)*Nl(1);
                         tempFq(1,2)*Nl(2) tempFq(2,2)*Nl(2);
                         tempFq(1,1)*Nl(2)+tempFq(1,2)*Nl(1) tempFq(2,1)*Nl(2)+tempFq(2,2)*Nl(1)];    
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
               J = dN(k)*Xe;
               me = me+rho*w(k)*fm(M(k)*Xe)*((N(k))'*N(k))*abs(det(J));
               ke = ke + (rho*fk(M(k)*Xe)*(cs/b_scal)^2*w(k)*((N(k))'*N(k))*abs(det(J)) +  w(k)*(Bt'*D./dt*Bq)*abs(det(J)));
               ce = ce + (rho*fc(M(k)*Xe)*(cs/b_scal)*w(k)*((N(k))'*N(k))*abs(det(J)) + w(k)*(Bt'*D./dt*Beps)*abs(det(J))); 
           end
   end   
   % Define global address vector ig( ) for element dofs.
   ig  = zeros(1,nne*dof);
   for ii = 1:nne
     for jj = 1:dof
       ig((ii-1)*dof+jj) = (Te(ii)-1)*dof + jj;
     end
   end
   % Add element contributions.
   for ii = 1:nne*dof
     for jj = 1:nne*dof
       Mass(ig(ii),ig(jj)) = Mass(ig(ii),ig(jj)) + me(ii,jj);
       Stiff(ig(ii),ig(jj)) = Stiff(ig(ii),ig(jj)) + ke(ii,jj);
       Damp(ig(ii),ig(jj)) = Damp(ig(ii),ig(jj)) + ce(ii,jj);
     end
   end
end



end

