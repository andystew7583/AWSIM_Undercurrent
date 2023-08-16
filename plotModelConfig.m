%%%
%%% plotModelConfig.m
%%%
%%% Plots the wind forcing, offshore stratification and bathymetry for our
%%% undercurrent test cases.
%%% 

%%% Load constant parameters 
constants;   

%%% Physical parameters
f0 = 8e-5;                   %%% Coriolis parameter 
beta = 0e-11;               %%% Coriolis parameter gradient
Ly = 400*m1km;                %%% Domain length.   
Lx = 2*Ly;
H = 4000;                     %%% Ocean depth  
Cd = 2e-3;                %%% Quadratic bottom drag
tau0 = 0.05;                  %%% Wind stress maximum      
tRelax_u = 1*t1day;         %%% Sponge relaxation time scales
tRelax_v = 1*t1day;
tRelax_h = 1*t1day;          
Lrelax = 50*m1km;           %%% Sponge relaxation width

%%% Stratification parameters
rho0 = 1000;                  %%% Reference density 
dz = 1;                       %%% Discretization step for continuous density profile
zz = 0:-dz:-H;                %%% Grid for continuous density profile
hs = 400;                     %%% Exponential decay scale for density profile
hl = 40000;                   %%% Linear decay scale for density profile 
rhoBot = 1028;                %%% Sea floor density in continuous density profile
rhoSurf = 1024;               %%% Sea surface density in continuous density profile
rhoRange = rhoBot-rhoSurf;

%%% Grids
Ny = 128;  
Nlay = 7;
d = Ly/Ny;  
Nx = round(Lx/d);
xx_q = 0:d:Lx;
yy_q = 0:d:Ly;  
xx_h = d/2:d:Lx-d/2;
yy_h = d/2:d:Ly-d/2;    
[YY_q,XX_q] = meshgrid(yy_q,xx_q);
[YY_h,XX_h] = meshgrid(yy_h,xx_h);
[YY_u,XX_u] = meshgrid(yy_h,xx_q(1:Nx));
[YY_v,XX_v] = meshgrid(yy_q(1:Ny),xx_h);    

%%% Will be incremented as more diagnostic figures are produced
fignum = 0;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% DENSITY STRATIFICATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Exponential/inear continuous density profile
rhoRef = rhoBot - rhoRange*(exp(zz/hs)+zz/hl-exp(-H/hs)+H/hl)/(1-exp(-H/hs)+H/hl); 

figure(1);
plot(rhoRef,zz);
xlabel('Offshore potential density profile (kg/m^3)');
ylabel('Elevation (m)');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% BOTTOM TOPOGRAPHY %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           

%%% Exponentials and linear slope with matched gradients
smax = 1.5e-1; %%% Max slope steepness
Ws_shelf = 5*m1km; %%% Reference slope exponential width for continental shelf
Ws_deep = 20*m1km; %%% Reference slope exponential width for deep ocean
Ys_shelf = 75*m1km; %%% Latitude of center of continental slope
H_shelf = 100; %%% Continental shelf water column thickness  
DY_canyons = 25*m1km; %%% Meridional amplitude of slope canyons
N_canyons = 4; %%% Number of slope canyons
X_canyons = Lx/4-Lx/N_canyons/2; %%% Longitude of first slope canyon
H_ridge = 0; %%% Meridional ridge height
W_ridge = 50*m1km; %%% Zonal width of meridional ridge

%%% Matrices defining spatially-dependent geometric parameters
Zbot = -H * ones(Nx,Ny); %%% Flat reference ocean floor depth    
Zshelf = -H_shelf * ones(Nx,Ny); %%% Continental shelf elevation
Ws_upper = Ws_shelf * ones(Nx,Ny); %%% Exponential width of upper slope
Ws_lower = Ws_deep * ones(Nx,Ny); %%% Exponential width of lower slope, reverts to Ws_shelf where ridge elevation is high
Delta_upper = Ws_upper * smax; %%% Exponential elevation change of upper slope
Delta_lower = Ws_lower * smax; %%% Exponential elevation change of lower slope
Zs_upper = (Zshelf - Delta_upper); %%% Elevation of upper/mid slope join
Zs_lower = (Zbot + Delta_lower); %%% Elevation of lower/mid slope join  
Ys_upper = Ys_shelf - DY_canyons*cos(2*pi*N_canyons*(XX_h-X_canyons)/Lx); %%% Latitude of upper/mid slope join
Ys_lower = Ys_upper + ((Zs_upper-Zs_lower)/smax); %%% Latitude of lower/mid slope join

%%% Construct the bathymetry by joining exponential upper and lower
%%% slopes with a linear mid-slope
etab = zeros(Nx,Ny);
idx_upper = YY_h < Ys_upper;
idx_mid = (YY_h >= Ys_upper) & (YY_h < Ys_lower);
idx_lower = YY_h > Ys_lower;    
etab_upper = Zshelf - Delta_upper .* exp((YY_h-Ys_upper)./Ws_upper);
etab_mid = Zs_upper - smax.*(YY_h-Ys_upper);
etab_lower = Zbot + Delta_lower .* exp(-(YY_h-Ys_lower)./Ws_lower);
etab(idx_upper) = etab_upper(idx_upper);
etab(idx_mid) = etab_mid(idx_mid);
etab(idx_lower) = etab_lower(idx_lower);

%%% Add random topographic variations
etab_lambda = 100*m1km;
etab_E0 = 1;
etab_amp = 200;
%   etab_amp = 500;
etab_pert = genRandIC (etab_lambda,etab_E0,Nx,Ny,Ly);
etab_pert = etab_amp*etab_pert/sqrt(mean(etab_pert(:).^2));
etab_pert = etab_pert.*abs(etab)/H;
%   etab_pert = etab_pert/2-min(min(etab_pert));
etab_pert = etab_pert.*(abs(etab)/H).^2;  
etab = etab + etab_pert;

%%% Plot topography
figure(2);
surf(XX_h/1000,YY_h/1000,etab);
shading interp;
colorbar;
pbaspect([2,1,1]);
xlabel('Alongshore distance (km)');
ylabel('Offshore distance (km)');
zlabel('Elevation (m)');

figure(3)
plot(yy_h/1000,etab(Nx/2,:));
hold on;
plot(yy_h/1000,etab(3*Nx/8,:));
hold off;  
xlabel('Offshore distance (km)');
zlabel('Elevation (m)');



%%%%%%%%%%%%%%%%%%%%%%%
%%%%% WIND STRESS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%

%%% Coast-localized wind
Lwind = 200*m1km;
taux = -tau0*sin(pi*YY_u/Lwind).^2 / rho0;  
taux(YY_u>Lwind) = 0;

%%% Realistic wind shape with coastal drop-off
%   Lwind = 75*m1km;
%   taux = tau0*tanh(YY_u/Lwind).^2 / rho0;

figure(4)
plot(yy_h/1000,taux(1,:)*rho0);
title('Wind stress (N/m^2)');
xlabel('Offshore distance (km)');
zlabel('Alongshore wind stress (N/m^2)');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% RELAXATION SETUP %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Relaxation time scale
relTime = -ones(Nx,Ny);    
relTime(YY_h>Ly-Lrelax) = 1 ./ (1-(Ly-YY_h(YY_h>Ly-Lrelax))/Lrelax);
eTime = tRelax_h * repmat(reshape(relTime,[1 Nx Ny]),[Nlay 1 1]);
uTime = tRelax_u * repmat(reshape(relTime,[1 Nx Ny]),[Nlay 1 1]);
vTime = tRelax_v * repmat(reshape(relTime,[1 Nx Ny]),[Nlay 1 1]);
hTime = tRelax_h * relTime;

figure(5);
plot(YY_h(1,:)/1000,squeeze(eTime(:,Nx/2,:))/t1day);
title('Relaxation time scales');
xlabel('Offshore distance (km)');
ylabel('Relaxation time scale (days)')
  