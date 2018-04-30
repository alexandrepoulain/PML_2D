function [outputArg1,outputArg2] = assemble_M_global(inputArg1,inputArg2)
% This function assembles the mass matrix for the 2D PML
% Inputs: rho: density of the PML
%         globX: coordinates of the nodes in PML
%         globT: connectivity in PML
%         a0: coefficient of attenuation for evanescent waves
%         x0: beginning of the PML
%         L2: length of the PML in x direction
%         h: length of the medium in y direction
%         n: order of the attenuation functions
%         dof: deggres of freedom for a node
% Function of attenuation in 2D
fe1 = @(x) a0(1)*((x-x0)/L2)^n; 
fe2 = @(y) a0(2)*((y-h)/L2)^n;
% function fm
fm = @(x) (1+fe1(x(1)))*(1+fe2(x(2)));
% number of nodes in an element
nn = size(globT,2)-1;
% Gauss quadrature 2 points (In each direction)
w = [1 1 1 1];
eta = [-1/sqrt(3) -1/sqrt(3) 1/sqrt(3) 1/sqrt(3)];
ksi = [-1/sqrt(3) 1/sqrt(3) 1/sqrt(3) -1/sqrt(3)];

% For each element in the PML (each row in globT)
MB = zeros(size(globX,1)*dof,size(globX,1)*dof);
for ii = 1:size(globT,1)
   % Get back the position of the nodes
   Xe = globX(globT(ii,1:nn),:);
   
   me = zeros(dof*nn,dof*nn);
   % 3 cases: Med, PML or interface
   
   % For each node
   for k=1:nn
       paramXe = paramX(globXe(1,k),:);
       switch paramXe{1}
           case 'MED' 
               dN = [ -(1-eta(k))  (1-eta(k))  (1+eta(k)) -(1+eta(k));
	           -(1-ksi(k)) -(1+ksi(k))  (1+ksi(k))  (1-ksi(k)) ]/4;
               local_N = N(ksi(k),eta(k));
               J = dN*Xe;
               me = me + rho*w(k)*( local_N'*local_N )*abs(det(J));
           case 'PML'
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
               me = me+rho*w(k)*fm(M*Xe)*(N'*N)*abs(det(J));
           case 'inter'
   end
   % assemble
   % Define global address vector ig( ) for element dofs.
   Te = globT(ii,:);
   ig  = zeros(1,nn*dof);
   for i = 1:nn
     for j = 1:dof
       ig((i-1)*dof+j) = (Te(i)-1)*dof + j;
     end
   end

   % Add element contributions.
   for i = 1:nn*dof
    for j = 1:nn*dof
       MB(ig(i),ig(j)) = MB(ig(i),ig(j)) + me(i,j);
     end
   end
end


end

