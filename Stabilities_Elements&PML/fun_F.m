function [Mat_FQ, Mat_Feps] = fun_F(a0,b0,x0,h,L2,n,cs,b_scal,dt,y,Xe)
% This function calculates the Feps amd FQ matrices
% Inputs: a0,b0 : coefficients of attenuation in both directions
%         x0: start of PML in x direction
%         y0: start of PML in y direction
%         L2: Length of PML (in x direction)
%         h:  Length of PML (in y direction)
%         n: order of attenuation on PML
%         cs: celerity of S-wave
%         b_scal: scalar
%         dt: time step on PML
%         y: local position of Gauss Point

% rotation matrix
theta = 0;
Q = [cos(theta) sin(theta); -sin(theta) cos(theta)];
% Function of attenuation in 2D
% fe1 = @(x) a0(1)*((x-x0)/L2)^n * (x>x0); 
% fe2 = @(y) a0(2)*((y-h)/L2)^n * (y>h);
% fp1 = @(x) b0(1)*((x-x0)/L2)^n*(x>x0); 
% fp2 = @(y) b0(2)*((y-h)/L2)^n*(y>h);
% stability
fe1 = @(x) a0(1); 
fe2 = @(y) a0(2);
fp1 = @(x) b0(1); 
fp2 = @(y) b0(2);
% Matrix Fp, Fe, Fte, Ftp
Fe = @(x) Q'*[1+fe1(x(1)) 0; 0 1+fe2(x(2))]*Q;
Fp = @(x) Q'*[(cs/b_scal)*fp1(x(1)) 0; 0 (cs/b_scal)*fp2(x(2))]*Q;
Fte = @(x) Q'*[1+fe2(x(2)) 0; 0 1+fe1(x(1))]*Q;
Ftp = @(x) Q'*[(cs/b_scal)*fp2(x(2)) 0; 0 (cs/b_scal)*fp1(x(1))]*Q;
% Fe = @(x) Q'*[1+fe1(x(1)) 0; 0 1]*Q;
% Fp = @(x) Q'*[(cs/b_scal)*fp1(x(1)) 0; 0 0]*Q;
% Fte = @(x) Q'*[1 0; 0 1+fe1(x(1))]*Q;
% Ftp = @(x) Q'*[0 0; 0 (cs/b_scal)*fp1(x(1))]*Q;
% Matrix Fl, Fq
Fl = @(x)(Fp(x)+(Fe(x)/dt))^(-1);
Feps = @(x) Fe(x)*Fl(x);
FQ = @(x) Fp(x)*Fl(x);
M = [((1-y(1))*(1-y(2)))  ((1+y(1))*(1-y(2)))   ((1+y(1))*(1+y(2)))  ((1-y(1))*(1+y(2))) ]/4;
pos = M*Xe;
tempFQ = FQ(pos); tempFeps = Feps(pos);
% Calculation of matrices
Mat_FQ= zeros(3,3);

Mat_FQ(1,1)=(tempFQ(1,1))^2;
Mat_FQ(1,2)=(tempFQ(2,1))^2;
Mat_FQ(1,3)=tempFQ(1,1)*tempFQ(2,1);
Mat_FQ(2,1)=(tempFQ(1,2))^2;
Mat_FQ(2,2)=(tempFQ(2,2))^2;
Mat_FQ(2,3)=tempFQ(1,2)*tempFQ(2,2);
Mat_FQ(3,1)=2*tempFQ(1,1)*tempFQ(1,2);
Mat_FQ(3,2)=2*tempFQ(2,1)*tempFQ(2,2);
Mat_FQ(3,3)=tempFQ(1,1)*tempFQ(2,2)+tempFQ(1,2)*tempFQ(2,1);

Mat_Feps= zeros(3,3);

Mat_Feps(1,1)=(tempFeps(1,1))^2;
Mat_Feps(1,2)=(tempFeps(2,1))^2;
Mat_Feps(1,3)=tempFeps(1,1)*tempFeps(2,1);
Mat_Feps(2,1)=(tempFeps(1,2))^2;
Mat_Feps(2,2)=(tempFeps(2,2))^2;
Mat_Feps(2,3)=tempFeps(1,2)*tempFeps(2,2);
Mat_Feps(3,1)=2*tempFeps(1,1)*tempFeps(1,2);
Mat_Feps(3,2)=2*tempFeps(2,1)*tempFeps(2,2);
Mat_Feps(3,3)=tempFeps(1,1)*tempFeps(2,2)+tempFeps(1,2)*tempFeps(2,1);

end

