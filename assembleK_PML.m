function [KB] = assembleK_PML(rho,cs,b,XB,TB,a0,b0,x0,L2,h,n,dof)
% This function assembles the stiffness matrix for the 2D PML
% Inputs: rho: density of the PML
%         cs: celerity of shear wave
%         b: scalar factor 
%         XB: coordinates of the nodes in PML
%         TB: connectivity in PML
%         a0: coefficient of attenuation for evanescent waves
%         b0: coefficient of attenuation for propagating waves
%         x0: beginning of the PML
%         L2: length of the PML in x direction
%         h: length of the medium in y direction
%         n: order of the attenuation functions
%         dof: deggres of freedom for a node
% Function of attenuation in 2D
fp1 = @(x) b0(1)*((x-x0)/L2)^n*(x>x0); 
fp2 = @(y) b0(2)*((y-h)/L2)^n*(y>h);
% stability
% fp1 = @(x) b0(1); 
% fp2 = @(y) b0(2);
% function fc
fk = @(x) fp1(x(1))*fp2(x(2));
% number of nodes in an element
nn = size(TB,2)-1;
% Gauss quadrature 2 points (In each direction)
w = [1 1 1 1];
eta = [-1/sqrt(3) -1/sqrt(3) 1/sqrt(3) 1/sqrt(3)];
ksi = [-1/sqrt(3) 1/sqrt(3) 1/sqrt(3) -1/sqrt(3)];
% For each element in the PML (each row in TB)
KB = zeros(size(XB,1)*dof,size(XB,1)*dof);
for ii = 1:size(TB,1)
   % Get back the position of the nodes
   Xe = XB(TB(ii,1:nn),:);
   ke = zeros(dof*nn,dof*nn);
   % For each gauss  point
   for k=1:4
       % parametric derivatives
       dN = [-(1-eta(k))  (1-eta(k))  (1+eta(k)) -(1+eta(k)); 
	          -(1-ksi(k)) -(1+ksi(k))  (1+ksi(k))  (1-ksi(k))]/4;
       % Shape functions
       N = (1/4)*[((1-ksi(k))*(1-eta(k))) 0. ((1+ksi(k))*(1-eta(k))) 0. ((1+ksi(k))*(1+eta(k))) 0. ((1-ksi(k))*(1+eta(k))) 0. ; ...
                   0. ((1-ksi(k))*(1-eta(k))) 0. ((1+ksi(k))*(1-eta(k))) 0. ((1+ksi(k))*(1+eta(k))) 0. ((1-ksi(k))*(1+eta(k))) ] ;
       % Transform matrix 
       M = (1/4)*[((1-ksi(k))*(1-eta(k))) ((1+ksi(k))*(1-eta(k))) ((1+ksi(k))*(1+eta(k))) ((1-ksi(k))*(1+eta(k)))] ;
       % Global coordinates
       J = dN*Xe;
       ke = ke + rho*fk(M*Xe)*(cs/b)^2*w(k)*(N'*N)*abs(det(J));
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
       KB(ig(i),ig(j)) = KB(ig(i),ig(j)) + ke(i,j);
     end
   end
end

