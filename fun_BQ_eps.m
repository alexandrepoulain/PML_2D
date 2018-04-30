function [Mat_Bq, Mat_Beps] = fun_BQ_eps(eta,ksi,k,a0,b0,h,x0,L2,n,cs,b_scal,dt,Xe)
% This function calculates the matrices Bq and Beps
% Inputs: a0,b0 : coefficients of attenuation in both directions
%         x0: start of PML in x direction
%         y0: start of PML in y direction
%         L2: Length of PML (in x direction)
%         h:  Length of PML (in y direction)
%         n: order of attenuation on PML
%         cs: celerity of S-wave
%         b_scal: scalar
%         dt: time step on PML
%         Xe: global positions of thee 4 nodes of element
%         eta,ksi: The gauss points
%         k: at which gauss point we are

% material derivative
dN = @(eta,ksi)1/4*[-(1-eta)  (1-eta)  (1+eta) -(1+eta); 
	                -(1-ksi) -(1+ksi)  (1+ksi)  (1-ksi)];
% rotation matrix
theta = 0;
Q = [cos(theta) sin(theta); -sin(theta) cos(theta)];
% Function of attenuation in 2D
fe1 = @(x) a0(1)*((x-x0)/L2)^n * (x>x0); 
fe2 = @(y) a0(2)*((y-h)/L2)^n * (y>h);
fp1 = @(x) b0(1)*((x-x0)/L2)^n*(x>x0); 
fp2 = @(y) b0(2)*((y-h)/L2)^n*(y>h);

% stability
% fe1 = @(x) a0(1); 
% fe2 = @(y) a0(2);
% fp1 = @(x) b0(1); 
% fp2 = @(y) b0(2);
% Matrix Fp, Fe, Fte, Ftp
Fe = @(x) Q'*[1+fe1(x(1)) 0; 0 1+fe2(x(2))]*Q;
Fp = @(x) Q'*[(cs/b_scal)*fp1(x(1)) 0; 0 (cs/b_scal)*fp2(x(2))]*Q;
Fte = @(x) Q'*[1+fe2(x(2)) 0; 0 1+fe1(x(1))]*Q;
Ftp = @(x) Q'*[(cs/b_scal)*fp2(x(2)) 0; 0 (cs/b_scal)*fp1(x(1))]*Q;
% Matrix Fl, Fq
Fl = @(x)(Fp(x)+(Fe(x)/dt))^(-1);
Feps = @(x) Fe(x)*Fl(x);
Fq = @(x) Fp(x)*Fl(x);
M = (1/4)*[((1-ksi(k))*(1-eta(k))) ((1+ksi(k))*(1-eta(k))) ((1+ksi(k))*(1+eta(k))) ((1-ksi(k))*(1+eta(k)))] ;
pos = M*Xe; % position of the node in natural coordinates 
% calculation of the matrices
tempFl = Fl(pos);
tempdN = dN(eta(k),ksi(k));
J = tempdN*Xe;
tempdN = J\tempdN;
tempFesp = Feps(pos);tempFq = Fq(pos);
% calculation of Nl,Bq,Nte and Ntp
ind = 1;
Mat_Beps = zeros(3,8);Mat_Bq = zeros(3,8);
for kk=1:4
    Nl = [tempFl(1,1)*tempdN(1,kk)+tempFl(1,2)*tempdN(2,kk);
         tempFl(2,1)*tempdN(1,kk)+tempFl(2,2)*tempdN(2,kk)];
    Mat_Beps(:,ind:ind+1) = [tempFesp(1,1)*Nl(1) tempFesp(2,1)*Nl(1);
         tempFesp(1,2)*Nl(2) tempFesp(2,2)*Nl(2);
         tempFesp(1,1)*Nl(2)+tempFesp(1,2)*Nl(1) tempFesp(2,1)*Nl(2)+tempFesp(2,2)*Nl(1)];
    % calculation of Nl,Bq,Nte and Ntp
    Nl = [tempFl(1,1)*tempdN(1,kk)+tempFl(1,2)*tempdN(2,kk);
         tempFl(2,1)*tempdN(1,kk)+tempFl(2,2)*tempdN(2,kk)];
    Mat_Bq(:,ind:ind+1) = [tempFq(1,1)*Nl(1) tempFq(2,1)*Nl(1);
         tempFq(1,2)*Nl(2) tempFq(2,2)*Nl(2);
         tempFq(1,1)*Nl(2)+tempFq(1,2)*Nl(1) tempFq(2,1)*Nl(2)+tempFq(2,2)*Nl(1)];
    ind=ind+2;
end
end

