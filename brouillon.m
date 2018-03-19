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