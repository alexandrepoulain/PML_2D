function [M] = assembleM(XA,TA,dof,Nef,rho)
% Subroutine to assemble the mass matrix
% Inputs: XA: coordinates of the nodes
%         TA: Connectivity matrix
%         dof: degree of freedom for a node
%         Nef: number of finite elements
%         rho: density
% Gauss Quadrature 2 points in each direction
xi1 = -1/sqrt(3);
xi2 = 1/sqrt(3);
eta1 = -1/sqrt(3);
eta2 = 1/sqrt(3);
w1 = 1;
w2 = 1;
% Nodal shape functions
N1 = @(xi, eta) 1/4*(1-xi)*(1-eta);
N2 = @(xi, eta) 1/4*(1+xi)*(1-eta);
N3 = @(xi, eta) 1/4*(1+xi)*(1+eta);
N4 = @(xi, eta) 1/4*(1-xi)*(1+eta);
N = @(xi,eta) [N1(xi,eta) 0 N2(xi,eta) 0 N3(xi,eta) 0 N4(xi, eta) 0;
    0 N1(xi,eta) 0 N2(xi,eta) 0 N3(xi, eta) 0 N4(xi,eta)];

% Strain matrix
B = @(xi,eta) [-(1-eta)/4 0 (1-eta)/4 0 (1+eta)/4 0 -(1+eta)/4 0 ;
    0 -(1-xi)/4 0 -(1+xi)/4 0 (1+xi)/4 0 (1-xi)/4 ;
    -(1-xi)/4 -(1-eta)/4 -(1+xi)/4 (1-eta)/4 (1+xi)/4 (1+eta)/4 (1-xi)/4 -(1+eta)/4];
% Derivative of N
dN = @(xi,eta)[-(1-xi)/4 -(1-eta)/4 -(1+xi)/4 (1-eta)/4;
    (1+xi)/4 (1+eta)/4 (1-xi)/4 -(1+eta)/4];


M = zeros(size(XA,1)*dof, size(XA,1)*dof);
% For each element
for i = 1:Nef
   % Position of the nodes
   Xe = XA(TA(i,1:4),:);
   Te = TA(i,:);
   % Calculation of submatrix
   me = rho*(w1*(w1*(N(xi1,eta1))'*N(xi1,eta1)*abs(det(dN(xi1,eta1)*Xe))+w2*(N(xi1,eta2))'*N(xi1,eta2))*abs(det(dN(xi1,eta2)*Xe))+...
       w2*(w1*(N(xi2,eta1))'*N(xi2,eta1)*abs(det(dN(xi2,eta1)*Xe))+w2*(N(xi2,eta2))'*N(xi2,eta2)*abs(det(dN(xi2,eta2)*Xe))));
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

