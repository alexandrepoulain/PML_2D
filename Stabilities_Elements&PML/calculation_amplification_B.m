function [B] = calculation_amplification_B(dt,beta,gamma,bulk,mu,ng,TB,...
    XB,L2,n,a0,b0,x0,h,cs,b_scal,M,C,Ceff,K,Keff,nnodes,dof)
% Function to calculate the amplification matrix A
% Inputs: dt: time step
%         beta, gamma: Newmark parameters
%         bulk,mu: bulk modulus and first lame coefficient
%         ng: order of quadrature
%         TB, XB: connectivity and coordinates matrices
%         L2, h: length of PML in x and y directions
%         n: order of attenuation
%         a0,b0: attenuation coefficients from kukucoban
%         x0: start of PML
%         cs: S-wave celerity
%         b_scal: scalar
%         M,C,Ceff,K,Keff: mass, damping, effective damping, rigidity,
%         effective rigidity matrices


% matrices to multiply it with eps, Eps, Sig, depending on the order of the
% quadrature
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
% constitutive matrix plane strain
D  = [bulk+4/3*mu bulk-2*mu/3 0
	  bulk-2*mu/3 bulk+4/3*mu 0
	    0      0        mu]; 
% material derivative
dN = @(k)1/4*[-(1-eta(k))  (1-eta(k))  (1+eta(k)) -(1+eta(k)); 
	          -(1-ksi(k)) -(1+ksi(k))  (1+ksi(k))  (1-ksi(k))];
% For each element
for ii = 1:size(TB,1)
    nn = size(TB,2)-1; % number of nodes
    Xe = XB(TB(ii,1:nn),:); % position of the 4 nodes
    % For each gauss poimt in the element
    for k = 1:ng^2
        y=[ksi(k) eta(k)];
        tempdN = dN(k);
        J = tempdN*Xe;
        Mat_D = D/dt;
        [Mat_B, Mat_Bp] = fun_B(a0,b0,L2,n,x0,h,cs,b_scal,Xe,eta,ksi,dt,k);
        [Mat_FQ, Mat_Feps] = fun_F(a0,b0,x0,h,L2,n,cs,b_scal,dt,y,Xe);
        [Mat_BQ, Mat_Beps] = fun_BQ_eps(eta,ksi,k,a0,b0,h,x0,L2,n,cs,b_scal,dt,Xe);

        Cell_epsi_n{k}((ii-1)*8+1:(ii-1)*8+8,(ii-1)*3+1:(ii-1)*3+3) = Mat_B' * Mat_D * Mat_Feps / dt * abs(det(J))*w(k);
        Cell_Epsi_n{k}((ii-1)*8+1:(ii-1)*8+8,(ii-1)*3+1:(ii-1)*3+3) = Mat_B' * Mat_D * Mat_FQ * abs(det(J))*w(k);
        Cell_Sig_n{k}((ii-1)*8+1:(ii-1)*8+8,(ii-1)*3+1:(ii-1)*3+3) = Mat_Bp' * abs(det(J))*w(k);
        Cell_Bepsi_n{k}(:,:,ii) = Mat_Beps / dt;
        Cell_BQ_n{k}(:,:,ii) = Mat_BQ  / dt;
        Cell_Fepsi_n{k}((ii-1)*3+1:(ii-1)*3+3,(ii-1)*3+1:(ii-1)*3+3) = Mat_Feps / (dt^2);
        Cell_FQ_n{k}((ii-1)*3+1:(ii-1)*3+3,(ii-1)*3+1:(ii-1)*3+3) = Mat_FQ / dt;
        Cell_Dsig{k}((ii-1)*3+1:(ii-1)*3+3,(ii-1)*3+1:(ii-1)*3+3) = Mat_D * dt;
    end
end
% reorganisation not into cells but matrices
for kk = 1:ng^2
    Mat_epsi_n(:,(kk-1)*3+1:(kk-1)*3+3) = Cell_epsi_n{kk}(:,:);
    Mat_Epsi_n(:,(kk-1)*3+1:(kk-1)*3+3) = Cell_Epsi_n{kk}(:,:);
    Mat_Sig_n(:,(kk-1)*3+1:(kk-1)*3+3)  = Cell_Sig_n{kk}(:,:);
    Mat_Bepsi_n((kk-1)*3+1:(kk-1)*3+3,:) = Cell_Bepsi_n{kk}(:,:);
    Mat_BQ_n((kk-1)*3+1:(kk-1)*3+3,:) = Cell_BQ_n{kk}(:,:);
    Mat_Fepsi_n((kk-1)*3+1:(kk-1)*3+3,(kk-1)*3+1:(kk-1)*3+3) = Cell_Fepsi_n{kk}(:,:);
    Mat_FQ_n((kk-1)*3+1:(kk-1)*3+3,(kk-1)*3+1:(kk-1)*3+3) = Cell_FQ_n{kk}(:,:);
    Mat_Dsig((kk-1)*3+1:(kk-1)*3+3,(kk-1)*3+1:(kk-1)*3+3) = Cell_Dsig{kk}(:,:);
end
% first line
B(1:nnodes*dof,:) = [-eye(size(M))+dt*(1-gamma)*M^-1*(C+Ceff) ...
    dt*(1-gamma)*M^-1*(K+Keff) ...
    dt*(1-gamma)*M^-1*Mat_epsi_n ...
    -dt*(1-gamma)*M^-1*Mat_Epsi_n ...
    dt*(1-gamma)*M^-1*Mat_Sig_n];
B(nnodes*dof+1:2*nnodes*dof,:) = [-dt*eye(size(M))+dt^2*(0.5-beta)*M^-1*(C+Ceff) ...
    -eye(size(M))+dt^2*(0.5-beta)*M^-1*(K+Keff) ...
    dt^2*(0.5-beta)*M^-1*Mat_epsi_n ...
    -dt^2*(0.5-beta)*M^-1*Mat_Epsi_n ...
    dt^2*(0.5-beta)*M^-1*Mat_Sig_n];
B(2*nnodes*dof+1:2*nnodes*dof+(3*ng^2),:) = [-Mat_Bepsi_n ...
    -Mat_BQ_n ...
    -Mat_Fepsi_n ...
    Mat_FQ_n ...
    zeros(size(Mat_Fepsi_n))];
B(2*nnodes*dof+(3*ng^2)+1:2*nnodes*dof+2*(3*ng^2),:) = [zeros(size(Mat_BQ_n)) ...
    zeros(size(Mat_BQ_n)) ...
    zeros(size(Mat_Fepsi_n)) ...
    -eye(size(Mat_Fepsi_n)) ...
    zeros(size(Mat_Fepsi_n))];
B(2*nnodes*dof+2*(3*ng^2)+1:2*nnodes*dof+3*(3*ng^2),:) = [zeros(size(Mat_BQ_n)) ...
    zeros(size(Mat_BQ_n)) ...
    zeros(size(Mat_Fepsi_n)) ...
    zeros(size(Mat_Fepsi_n)) ...
    -eye(size(Mat_Fepsi_n))];
end

