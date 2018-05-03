%% PML 2D enhanced (2nd order) Lamb's Test
% Parameters of input Ricker Wave
tp = 3; % Fundamental period
ts = 3; % Time shift
A = 1; % Amplitude

% Material parameters
rho_1 = 1; % density
mu_1 = 1; % 1st Lamé Coefficient 
nu_1 = 0.25; % Poisson modulus

cs = (mu_1/rho_1)^0.5 ; % Shear wave velocity
lmbd_s = cs*tp; % wave length of shear wave

L = 4/3*lmbd_s; % Size of the soil
