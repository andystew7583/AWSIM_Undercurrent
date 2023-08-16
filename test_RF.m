%%%
%%% test_RF.m
%%%
%%% Tests random forcing formulation.
%%%

%%% Load constant parameters 
constants;   

%%% Physical parameters
f0 = 8e-5;                   %%% Coriolis parameter 
beta = 0e-11;               %%% Coriolis parameter gradient
Ly = 400*m1km;                %%% Domain length.   
Lx = 2*Ly;
Nx = 512;
Ny = 256;
dx = Lx/Nx;
dy = Ly/Ny;
H = 4000;                     %%% Ocean depth  
Lrelax = 50*m1km;

%%% Random forcing parameters
useBaroclinicForcing = true; %%% Apply random forcing to near-surface layers only  
RF_tau = 864000; %%% Autocorrelation timescale  
lambdaK = 80*m1km; %%% Peak energy input lengthscale for random forcing
Ymin_RF = 0; %%% No random forcing below this latitude
Ymax_RF = Ly-Lrelax; %%% No random forcing above this latitude
Lramp_RF = 30*m1km; %%% Latitudinal width of region over which to ramp up the random forcing
  
%%% Generate forcing mask in spectral space
k = [0:1:Nx/2-1,-Nx/2:1:-1]; %%% Zonal wavenumber
K_xk = 2*pi.*(k)./Lx;
l = [0:1:Ny/2-1,-Ny/2:1:-1]; %%% Meridional wavenumber
K_yl = 2*pi.*(l)./Ly;
[K_ykl,K_xkl] = meshgrid(K_yl, K_xk);
K = sqrt(K_xkl.^2 + K_ykl.^2); %%% Absolute wavenumber 
K_0 = 2*pi/lambdaK; %%% Most strongly forced wavenumber
W = K_0/8; %%% Exponential width of energy band in wavenumber space   
RF_mask_fft = 2*exp(-((K-K_0)/W).^2);
% RF_mask_fft = ones(Nx,Ny);
RF_mask_norm = sqrt(0.5*sum(RF_mask_fft(:).^2));

RF_F0_rot = 0.75;
RF_kk_recip = 1.0 ./ K;
RF_kk_recip(1,1) = 0;
RFrot_fft = complex(zeros(Nx,Ny));
RFamp_rot_fft = RF_F0_rot .* RF_kk_recip .* RF_mask_fft ./ RF_mask_norm;

%%% Evolve over many forcing decorrelation time scales using the actual
%%% time step for this run to get a representative forcing field
%%% N.B. This code is basically taken directly
deltaT = 65;
Tratio = deltaT/RF_tau;
for n=1:20000
  RF_phase = 2*pi*rand(Nx,Ny);
  RFrot_fft = RFrot_fft - Tratio.*RFrot_fft + sqrt(2*Tratio).*exp(1i*RF_phase).*RFamp_rot_fft;
end
RFrot = Nx*Ny*real(ifft2(RFrot_fft));

%%% Test that RMS forcing magnitude matches specified forcing amplitude
dRFrot_dx = (RFrot(1:Nx,1:Ny)-RFrot([Nx 1:Nx-1],1:Ny))/dx;
dRFrot_dy = (RFrot(1:Nx,1:Ny)-RFrot(1:Nx,[Ny 1:Ny-1]))/dy;
gradsq_dRFrot = dRFrot_dx.^2+dRFrot_dy.^2;
sqrt(mean(gradsq_dRFrot(:)))