function [M] = assembleM(XA,TA,dof,Nef,rho)
% Subroutine to assemble the mass matrix
% Inputs: XA: coordinates of the nodes
%         TA: Connectivity matrix
%         dof: degree of freedom for a node
%         Nef: number of finite elements
%         rho: density
% Gauss Quadrature 2 points in each direction
w = [1 1 1 1];
eta = [-1/sqrt(3) -1/sqrt(3) 1/sqrt(3) 1/sqrt(3)];
ksi = [-1/sqrt(3) 1/sqrt(3) 1/sqrt(3) -1/sqrt(3)];
% Nodal shape functions
N1 = @(xi, eta) 1/4*(1-xi)*(1-eta);
N2 = @(xi, eta) 1/4*(1+xi)*(1-eta);
N3 = @(xi, eta) 1/4*(1+xi)*(1+eta);
N4 = @(xi, eta) 1/4*(1-xi)*(1+eta);
N = @(xi,eta) [N1(xi,eta) 0 N2(xi,eta) 0 N3(xi,eta) 0 N4(xi, eta) 0;
    0 N1(xi,eta) 0 N2(xi,eta) 0 N3(xi, eta) 0 N4(xi,eta)];

M = zeros(size(XA,1)*dof, size(XA,1)*dof);

% For each element
for i = 1:Nef
   % Position of the nodes
   Xe = XA(TA(i,1:4),:);
   nnodes = size(Xe,1);
   Te = TA(i,:);
   % Calculation of submatrix
   me = zeros(2*nnodes);
   for k=1:4
       dN = [ -(1-eta(k))  (1-eta(k))  (1+eta(k)) -(1+eta(k));
	           -(1-ksi(k)) -(1+ksi(k))  (1+ksi(k))  (1-ksi(k)) ]/4;
       local_N = N(ksi(k),eta(k));
       J = dN*Xe;
       me = me + rho*w(k)*( local_N'*local_N )*abs(det(J));
    end
   % Define global address vector ig( ) for element dofs.
   ne=4; dof=2;
   ig  = zeros(1,ne*dof);
   for ii = 1:ne
     for jj = 1:dof
       ig((ii-1)*dof+jj) = (Te(ii)-1)*dof + jj;
     end
   end
   % Add element contributions.
   for ii = 1:ne*dof
     for jj = 1:ne*dof
       M(ig(ii),ig(jj)) = M(ig(ii),ig(jj)) + me(ii,jj);
     end
   end
end

end

