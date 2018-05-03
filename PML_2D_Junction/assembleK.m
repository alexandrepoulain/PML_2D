function [K] = assembleK(XA,TA,dof,E,nu,Nef)
% Subroutine to assemble the stiffness matrix
% Inputs: XA : coordinates of the nodes
%         TA : connectivity
%         dof : degrees of freedom
%         E : Young modulus
%         nu : Poisson ratio
%         Nef : Number of elements
% Gauss Quadrature 2 points in each direction
w = [1 1 1 1];
eta = [-1/sqrt(3) -1/sqrt(3) 1/sqrt(3) 1/sqrt(3)];
ksi = [-1/sqrt(3) 1/sqrt(3) 1/sqrt(3) -1/sqrt(3)];
% initialization
K = zeros(size(XA,1)*dof, size(XA,1)*dof);
% Strain matrix
B = @(xi,eta) [-(1-eta)/4 0 (1-eta)/4 0 (1+eta)/4 0 -(1+eta)/4 0 ;
    0 -(1-xi)/4 0 -(1+xi)/4 0 (1+xi)/4 0 (1-xi)/4 ;
    -(1-xi)/4 -(1-eta)/4 -(1+xi)/4 (1-eta)/4 (1+xi)/4 (1+eta)/4 (1-xi)/4 -(1+eta)/4];

% constitutive matrix
D  = E/((1+nu)*(1-2*nu)) ...   
	*[ 1-nu   nu       0
	    nu   1-nu      0
	     0     0  (1-2*nu)/2 ];
     
for i = 1:Nef
   % Position of the nodes
   Xe = XA(TA(i,1:4),:);
   Te = TA(i,:);
   % for each gauss point
   nnodes = size(Xe,1);
   ke = zeros(2*nnodes);
   for k=1:4
       dN = [ -(1-eta(k))  (1-eta(k))  (1+eta(k)) -(1+eta(k)) 
	           -(1-ksi(k)) -(1+ksi(k))  (1+ksi(k))  (1-ksi(k)) ]/4;
       J = dN*Xe;  % Jacobian
       global_dN = J\dN; %Derivative of shape function in global coordinates
       % Construction of the strain matrix in global coordinates
       B  = [  global_dN(1,1)    0    global_dN(1,2)   0     global_dN(1,3)    0    global_dN(1,4)   0
	               0    global_dN(2,1)    0    global_dN(2,2)    0    global_dN(2,3)    0    global_dN(2,4)
	            global_dN(2,1) global_dN(1,1) global_dN(2,2) global_dN(1,2) global_dN(2,3) global_dN(1,3) global_dN(2,4) global_dN(1,4) ];
       ke = ke + w(k)*( B'*D*B )*abs(det(J));
   end
   % Define global address vector ig( ) for element dofs.
   ne=4;
   ig  = zeros(1,ne*dof);
   for ii = 1:ne
     for jj = 1:dof
       ig((ii-1)*dof+jj) = (Te(ii)-1)*dof + jj;
     end
   end
   % Add element contributions.
   for ii = 1:ne*dof
     for jj = 1:ne*dof
       K(ig(ii),ig(jj)) = K(ig(ii),ig(jj)) + ke(ii,jj);
     end
   end
end
end

