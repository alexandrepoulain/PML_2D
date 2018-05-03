function [Mat_epsi_n,Mat_Sig_n,vect_epsi_n,vect_epsi_n_2,vect_sig_n,...
    vect_Epsi_n,vect_Sig_n,vect_past,Mat_Bepsi_n,Mat_BQ_n,Mat_Fepsi_n,...
    Mat_FQ_n,Mat_Dsig,Mat_Epsi_n] = calculate_cste_matrices(Nef_2, ng,bulk,mu,TB,XB,...
    a0,b0,L2,n,x0,h,cs,b_scal,dt)
% This function calculates constant matrices used in 
% the calculation of stress, strain.  
% Inputs :  



% Initialization of sparse matrices (for performance)
for k=1:ng^2
    % In order to calculate the internal forces
    Mat_epsi_n{k} = sparse(Nef_2*8,Nef_2*3); % To calculate sigma multiplied by Bt
    Mat_Sig_n{k} = sparse(Nef_2*8,Nef_2*3); % same
    vect_epsi_n{k} = zeros(Nef_2*3,1);  % epsilon at each gauss point
    vect_epsi_n_2{k} = zeros(Nef_2*3,1); % first part of the calculation of epsilon
    vect_sig_n{k} = zeros(Nef_2*3,1); % stress at each gauss point 
    vect_Epsi_n{k} = zeros(Nef_2*3,1); % integral of the strain at each gauss point 
    vect_Sig_n{k} = zeros(Nef_2*3,1); % integral of the stress at each gauss point
    vect_past = zeros(Nef_2*8,1); % Vector of internal forces
    Mat_Bepsi_n{k} = zeros(3,8,Nef_2); % To calculate Besp/dt (in the calculation of the strain at the next time step)
    Mat_BQ_n{k} = zeros(3,8,Nef_2); % Bq/dt same
    Mat_Fepsi_n{k} = sparse(Nef_2*3,Nef_2*3); % To calculate Feps/dt²
    Mat_FQ_n{k} = sparse(Nef_2*3,Nef_2*3); % To calculate Fq/dt
    Mat_Dsig{k} = sparse(Nef_2*3,Nef_2*3); % To calculate D/dt (in the calculation of the stress at the next time step)
end
% Depending on the order of the quadrature used we need different Gauss
% points and weights
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
% material derivative
dN = @(k)1/4*[-(1-eta(k))  (1-eta(k))  (1+eta(k)) -(1+eta(k)); 
	                -(1-ksi(k)) -(1+ksi(k))  (1+ksi(k))  (1-ksi(k))];
% constitutive matrix plane strain
D  = [bulk+4/3*mu bulk-2*mu/3 0
	  bulk-2*mu/3 bulk+4/3*mu 0
	    0      0        mu]; 
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
        Mat_epsi_n{k}((ii-1)*8+1:(ii-1)*8+8,(ii-1)*3+1:(ii-1)*3+3) = Mat_B' * Mat_D * Mat_Feps / dt * abs(det(J))*w(k);
        Mat_Epsi_n{k}((ii-1)*8+1:(ii-1)*8+8,(ii-1)*3+1:(ii-1)*3+3) = Mat_B' * Mat_D * Mat_FQ * abs(det(J))*w(k);
        Mat_Sig_n{k}((ii-1)*8+1:(ii-1)*8+8,(ii-1)*3+1:(ii-1)*3+3) = Mat_Bp' * abs(det(J))*w(k);
        Mat_Bepsi_n{k}(:,:,ii) = Mat_Beps / dt;
        Mat_BQ_n{k}(:,:,ii) = Mat_BQ  / dt;
        Mat_Fepsi_n{k}((ii-1)*3+1:(ii-1)*3+3,(ii-1)*3+1:(ii-1)*3+3) = Mat_Feps / (dt^2);
        Mat_FQ_n{k}((ii-1)*3+1:(ii-1)*3+3,(ii-1)*3+1:(ii-1)*3+3) = Mat_FQ / dt;
        Mat_Dsig{k}((ii-1)*3+1:(ii-1)*3+3,(ii-1)*3+1:(ii-1)*3+3) = Mat_D * dt;
    end
end

end

