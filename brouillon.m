% Square elements
a = Lef;
b = Lef;
% Nodal shape functions
N1 = @(xi, eta) 1/4*(1-xi)*(1-eta);
N2 = @(xi, eta) 1/4*(1+xi)*(1-eta);
N3 = @(xi, eta) 1/4*(1+xi)*(1+eta);
N4 = @(xi, eta) 1/4*(1-xi)*(1+eta);

C = E_1/(1-nu_1) * [1 nu_1 0;
                    nu_1 1 0;
                    0 0 (1-nu_1)/2];

N = @(xi,eta) [N1(xi,eta) 0 N2(xi,eta) 0 N3(xi,eta) 0 N4(xi, eta) 0;
    0 N1(xi,eta) 0 N2(xi,eta) 0 N3(xi, eta) 0 N4(xi,eta)];

% Strain matrix
B = @(xi,eta) [-(1-eta)/4 0 (1-eta)/4 0 (1+eta)/4 0 -(1+eta)/4 0 ;
    0 -(1-xi)/4 0 -(1+xi)/4 0 (1+xi)/4 0 (1-xi)/4 ;
    -(1-xi)/4 -(1-eta)/4 -(1+xi)/4 (1-eta)/4 (1+xi)/4 (1+eta)/4 (1-xi)/4 -(1+eta)/4];
% Derivative of N
dN = @(xi,eta)[-(1-xi)/4 -(1-eta)/4 -(1+xi)/4 (1-eta)/4;
    (1+xi)/4 (1+eta)/4 (1-xi)/4 -(1+eta)/4];

%% Gauss Quadrature 2 points in each direction
xi1 = -1/sqrt(3);
xi2 = 1/sqrt(3);
eta1 = -1/sqrt(3);
eta2 = 1/sqrt(3);
w1 = 1;
w2 = 1;


%% Nodes with load
load_nodes = zeros(1,2);
compt = 1;
for i=1:nnA
  if (XA(i,1) == 0)
      load_nodes(compt) = i;
      compt = compt+1;
  end
end
Forn1  = zeros(dof*nnA,1) ; 
for i=1:length(load_nodes)
  Forn1((load_nodes(i)-1)*dof+1)=1 ;
end
%% Boundary conditions
% On SDA: fixed end
BCAx = [51 102];
BCAy = [51 102];



% for each element
for i=1:size(TB,1)
    nn = size(TB,2)-1; % number of total nodes
    Xe = XB(TB(i,1:nn),:); % coordinates of the nodes in the current element
    % for each gauss point
    for k=1:(ng^2)
        y=[ksi(k) eta(k)]; % position of the gauss point in natural coordinates 
        Jt = dN(eta(k),ksi(k))*Xe; % calculation of the jacobian of the mapping
        %% calculate B, Bp, Be for the current gauss point
        tempdN = [ -(1-y(2))  (1-y(2))  (1+y(2)) -(1+y(2)) 
	           -(1-y(1)) -(1+y(1))  (1+y(1))  (1-y(1)) ]/4;
        tempJt = tempdN*Xe ;
        tempdN = tempJt\tempdN ;
        temp_Fe= Fe(y);
        Fp=Fp(y);
        Be= zeros (3,8);
        Bp= zeros (3,8);
        for i=1:4
            Be(1,2*(i-1)+1)=temp_Fe(1,1)*dN(1,i)+temp_Fe(1,2)*dN(2,i);
            Be(2,2*(i-1)+2)=Fe(2,1)*dN(1,i)+Fe(2,2)*dN(2,i);
            Be(3,2*(i-1)+1)=Be(2,2*(i-1)+2);
            Be(3,2*(i-1)+2)=Be(1,2*(i-1)+1);
        end
        
        for i=1:4
            Bp(1,2*(i-1)+1)=Fp(1,1)*dN(1,i)+Fp(1,2)*dN(2,i);
            Bp(2,2*(i-1)+2)=Fp(2,1)*dN(1,i)+Fp(2,2)*dN(2,i);
            Bp(3,2*(i-1)+1)=Bp(2,2*(i-1)+2);
            Bp(3,2*(i-1)+2)=Bp(1,2*(i-1)+1);
        end            
        B=Be+dt*Bp;
        %% Matrices Bq end Beps
        Beps{k}(:,ind:ind+1,i) = [tempFesp(1,1)*Nl(1) tempFesp(2,1)*Nl(1);
             tempFesp(1,2)*Nl(2) tempFesp(2,2)*Nl(2);
             tempFesp(1,1)*Nl(2)+tempFesp(1,2)*Nl(1) tempFesp(2,1)*Nl(2)+tempFesp(2,2)*Nl(1)];
    end