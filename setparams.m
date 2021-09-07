%%%
%%% setparams.m
%%%
%%% Configures a shelf/slope channel with wind and buoyancy forcing.
%%% 
function setparams (local_home_dir,run_name)

  %%% Choose configuration
  %%% This affects the forcing that is permitted and the discretization in
  %%% the vertical.
  config = 'wind';
%   config = 'rand';
  
  %%% Set true to take averaged diagnostics
  full_diags = true;

  %%% Load constant parameters 
  constants;   
  
  %%% Run directory
  run_name = strtrim(run_name); 
  local_home_dir = strtrim(local_home_dir); 
  local_run_dir = fullfile(local_home_dir,run_name);
  mkdir(local_run_dir);
  pfname = fullfile(local_run_dir,[run_name,'_in']);   
  model_code_dir = fullfile('../AWSIM/',model_code_dir_name);
  
  %%% Cluster config
  %%% NOTE: You will need to edit matlab_common/createRunScript to add a
  %%% configuration for your cluster!
  use_cluster = false;
  use_intel = false;
  use_pbs = use_cluster;
  uname = 'astewart';
  cluster_addr = 'caolila.atmos.ucla.edu';
  cluster_home_dir = '/data2/astewart/AWSIM_WindAABW/runs';
  
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%% PARAMETERS %%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %%% To store parameters
  paramTypes;
  PARAMS = {};
  
  %%% Grid resolution
  N = 128;
  switch (config)
    case 'wind'
      Nlay = 7;
%       Nlay = 3;
%       Nlay = 4;
    case 'rand'
      Nlay = 4;
%       Nlay = 3;
  end
  
  %%% Physical parameters
  f0 = 8e-5;                   %%% Coriolis parameter 
  phiP = 0;                %%% Angle of rotation relative to north (0=N, pi/2=E, pi=S, 3pi/2=W), default 0
  beta = 0e-11;               %%% Coriolis parameter gradient
  beta_x = -beta*sin(phiP);     %%% Coriolis parameter gradient in x-direction
  beta_y = beta*cos(phiP);      %%% Coriolis parameter gradient in y-direction
  Ly = 400*m1km;                %%% Domain length.   
  Lx = 2*Ly;
  H = 4000;                     %%% Ocean depth  
  rb = 0;                    %%% Linear bottom drag
  Cd = 2e-3;                %%% Quadratic bottom drag
  Cd_surf = 0e-3;           %%% Surface drag coefficient
  switch (config)
    case 'wind'
      tau0 = 0.05;                  %%% Wind stress maximum    
      Fbaro0 = 0; %0.01;                 %%% Maximum depth-integrated along-slope pressure gradient force
    case 'rand'
      tau0 = 0;                  %%% Wind stress maximum    
      Fbaro0 = 0;                 %%% Maximum depth-integrated along-slope pressure gradient force
  end
  go_west = true;              %%% Set true for westward flow - reverses wind stress and initial flow    
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
  
  %%% Random forcing parameters
  switch (config)
    case 'wind'
      useRandomForcing = false;
      RF_F0_rot = 0; %%% Random forcing amplitude (m^2/s^2)
    case 'rand'
      useRandomForcing = true;
      RF_F0_rot = 0.75/rho0; %%% Random forcing amplitude (m^2/s^2)
  end  
  useBaroclinicForcing = true; %%% Apply random forcing to near-surface layers only  
  RF_tau = 10*t1day; %%% Autocorrelation timescale  
  lambdaK = 80*m1km; %%% Peak energy input lengthscale for random forcing
  Ymin_RF = 0; %%% No random forcing below this latitude
  Ymax_RF = Ly-Lrelax; %%% No random forcing above this latitude
  Lramp_RF = 30*m1km; %%% Latitudinal width of region over which to ramp up the random forcing
  
  %%% Temporal parameters  
  if (full_diags)
    tmax = 10*t1year;
  else
    tmax = 5*t1year;
  end
  savefreq = 5*t1day;   
  savefreqEZ = 0.1*t1day;  
  if (full_diags)
    savefreqAvg = t1year;
    savefreqUMom = t1year;
    savefreqVMom= t1year;
    savefreqThic = -1;
  else
    savefreqAvg = -t1year;
    savefreqUMom = -t1year;
    savefreqVMom= -t1year;
    savefreqThic = -1;
  end
  restart = 0;
  startIdx = 0;
  
  %%% Rigid lid-related parameters
  useRL = 1; %%% Set to 1 to use rigid lid, or 0 not to  
  use_MG = 1;
  tol = 1e-7;
  SOR_rp_max = 2.0;
  SOR_rp_min = 1.7;
  SOR_opt_freq = 1000; 
  
  %%% Grids
  Ny = N;  
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

  
  
  %%% Option 1: uniform density jumps
  
%   %%% Offshore layer elevations
%   E0 = zeros(1,Nlay+1); 
%   E0(Nlay+1) = -H;
%   for k=2:Nlay
%     E0(k) = zz(find(rhoRef-rhoSurf>(k-1)*rhoRange/Nlay,1,'first')); %%% Break into Nlay-1 density jumps of equal sizes
%   end  
%   
%   %%% Offshore layer mid-depths and thicknesses
%   Zmid0 = 0.5*(E0(1:Nlay)+E0(2:Nlay+1));
%   H0 = E0(1:Nlay)-E0(2:Nlay+1) 
%   E0
%   
%   %%% Reduced gravities at layer interfaces   
%   gprime = g*rhoRange/(Nlay-1)/rho0; 
%   gg = [g gprime*ones(1,Nlay-1)]; 
  


  %%% Option 2: uniform layer thicknesses
  
%   %%% Offshore layer thicknesses, interface elevations, and mid-depths
%   H0 = H/Nlay*ones(1,Nlay);
%   E0 = zeros(1,Nlay+1); 
%   E0(2:Nlay+1) = -cumsum(H0);
%   Zmid0 = 0.5*(E0(1:Nlay)+E0(2:Nlay+1));
%   
%   rhoMean = zeros(1,Nlay);
%   for k=1:Nlay 
%     idx = find((zz<E0(k))  & (zz>E0(k+1)));
%     rhoMean(k) = mean(rhoRef(idx));
%   end
%   gg = [g g*diff(rhoMean)/rho0];
%   diff(rhoMean)
  
  



 

  %%% Option 3: Shallow upper layer then fixed spacing    
  switch (config)
  
    case 'wind'
    
      %%% Fixed spacing down from surface
      Hupper = 150; %%% For 7-layer case
      Hinner = 150;
%       Hupper = 250; %%% For 3-layer case
%       Hinner = 400;

    case 'rand' 
      
      Hupper = 50;
      Hinner = 200; %%% For 4-layer case
%       Hinner = 400; %%% For 3-layer case
  
  end    
   
  E0 = 0;
  if (Nlay >= 1)
    E0 = [E0 -Hupper];
  end
  if (Nlay >= 2)
    E0 = [E0 -Hupper-Hinner*(1:Nlay-2)];
  end  
  E0 = [E0 (-H)];
  H0 = E0(1:Nlay)-E0(2:Nlay+1);

  rhoMean = zeros(1,Nlay);
  for k=1:Nlay
    idx = find((zz<E0(k))  & (zz>E0(k+1)));
    rhoMean(k) = mean(rhoRef(idx));
  end
  gg = [g g*diff(rhoMean)/rho0];
  
  
  



%   %%% Option 4: Manual
%   
%   E0 = [0 -125 -250 -550 -H]; %%% For 4-layer case
%   H0 = E0(1:Nlay)-E0(2:Nlay+1);
%   rhoMean = zeros(1,Nlay);
%   for k=1:Nlay
%     idx = find((zz<E0(k))  & (zz>E0(k+1)));
%     rhoMean(k) = mean(rhoRef(idx));
%   end
%   gg = [g g*diff(rhoMean)/rho0];
  
  
  
  
  
  
  %%% Show stratification
  disp(H0)
  disp(diff(rhoMean))




  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% BOTTOM TOPOGRAPHY %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
%   %%% Put canyons in continental shelf
%   H_shelf = 100; %%% Continental shelf water column thickness  
%   W_canyons = 60*m1km;
%   H_canyons = 300;
%   X_canyons = [Lx/4 3*Lx/4];  
%   Z_shelf = -H_shelf * ones(Nx,Ny); %%% Continental shelf elevation
%   for m=1:length(X_canyons)
%     idx = (XX_h > X_canyons(m)-W_canyons/2) & (XX_h < X_canyons(m)+W_canyons/2);
%     Z_shelf(idx) = Z_shelf(idx) - H_canyons*0.5*(1+cos(2*pi*(XX_h(idx)-X_canyons(m))/(W_canyons))); 
%   end
%   
%   %%% Create continental slope
% %   Y_slope = 100*m1km; %%% Latitude of center of continental slope
% %   W_slope = 20*m1km;
% %   gam_h = 0.004;
% %   Y_scaled = (YY_h-Y_slope)/W_slope;     
% %   Z_slope = 0.5*(Z_shelf-H);
% %   H_slope = Z_shelf + H;
% %   etab = Z_slope + H_slope .* (0.25*sqrt((1-Y_scaled).^2 + 4*gam_h*Y_scaled.^2)-0.25*sqrt((1+Y_scaled).^2 +  4*gam_h*Y_scaled.^2)) / (1+4*gam_h)^(-1/2);
%   
%   %%% Exponentials and linear slope with matched gradients
%   smax = 1.5e-1; %%% Max slope steepness
%   Ws_shelf = 2*m1km; %%% Reference slope exponential width for continental shelf
%   Ws_deep = 12.5*m1km; %%% Reference slope exponential width for deep ocean
%   Ys_shelf = 80*m1km; %%% Latitude of center of continental slope  
%   
%   %%% Matrices defining spatially-dependent geometric parameters
%   Z_bot = -H * ones(Nx,Ny); %%% Flat reference ocean floor depth    
%   Ws_upper = Ws_shelf * ones(Nx,Ny); %%% Exponential width of upper slope
%   Ws_lower = Ws_deep * ones(Nx,Ny); %%% Exponential width of lower slope, reverts to Ws_shelf where ridge elevation is high
%   Delta_upper = Ws_upper * smax; %%% Exponential elevation change of upper slope
%   Delta_lower = Ws_lower * smax; %%% Exponential elevation change of lower slope
%   Zs_upper = (Z_shelf - Delta_upper); %%% Elevation of upper/mid slope join
%   Zs_lower = (Z_bot + Delta_lower); %%% Elevation of lower/mid slope join  
%   Ys_upper = Ys_shelf; %%% Latitude of upper/mid slope join
%   Ys_lower = Ys_upper + ((Zs_upper-Zs_lower)/smax); %%% Latitude of lower/mid slope join
%   
%   %%% Construct the bathymetry by joining exponential upper and lower
%   %%% slopes with a linear mid-slope
%   etab = zeros(Nx,Ny);
%   idx_upper = YY_h < Ys_upper;
%   idx_mid = (YY_h >= Ys_upper) & (YY_h < Ys_lower);
%   idx_lower = YY_h > Ys_lower;    
%   etab_upper = Z_shelf - Delta_upper .* exp((YY_h-Ys_upper)./Ws_upper);
%   etab_mid = Zs_upper - smax.*(YY_h-Ys_upper);
%   etab_lower = Z_bot + Delta_lower .* exp(-(YY_h-Ys_lower)./Ws_lower);a
%   etab(idx_upper) = etab_upper(idx_upper);
%   etab(idx_mid) = etab_mid(idx_mid);
%   etab(idx_lower) = etab_lower(idx_lower);
       

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

  %%% Uncomment to print useful geometric scalars
  Delta_upper(1,1)
  Delta_lower(1,1)
  Zs_upper(1,1)
  Zs_lower(1,1)
  Ys_upper(1,1)
  Ys_lower(1,1)
  
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
  
  %%% Add ridge
  if (H_ridge > 0)
    etab_ridge = H_ridge * exp(-((XX_h-Lx/2)/W_ridge).^2);
    etab_ridge = etab_ridge .* max(1-(etab+H)/H_ridge,0);
    etab = etab + etab_ridge;
  end
  
  %%% Plot topography
  fignum = fignum+1;
  figure(fignum);
  surf(XX_h,YY_h,etab);
  shading interp;
  colorbar;
%   plot(yy_h,etab(1,:))
  pbaspect([2,1,1]);

  fignum = fignum+1;
  figure(fignum);
  plot(yy_h,etab(Nx/2,:));
  hold on;
  plot(yy_h,etab(Nx/4,:));
  hold off;
  
  fignum = fignum+1;
  figure(fignum);
  plot(yy_q(2:end-1),diff(etab(Nx/2,:))/d);
  hold on;
  plot(yy_q(2:end-1),diff(etab(Nx/4,:))/d);
  hold off;
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% OTHER PARAMETERS %%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
   
  %%% Set timestep based on full gravity wave speed, assuming a that the
  %%% absolute velocity never exceeds the gravity wave speed. Then correct
  %%% dt so that tmax is an integer number of time steps.  
  c = calcWaveSpeed(H0,gg,useRL)
  Umax = 3;  
  switch (config)
    case 'wind'
      dt = 0.25*d/(c+Umax)  %%% Default time step choice
    case 'rand'
      dt = 0.25*d/(c+Umax)  %%% Default time step choice
%       dt = 0.2*d/(c+Umax)  %%% Default time step choice
  end  
%   dt = 0.04*d/(c+Umax) %%% For thinner Salmon layers (h0=1)
%   dt = 0.01*d/(c+Umax) %%% For even thinner Salmon layers (h0=0.5)
%   dt = 0.001*d/(c+Umax) %%% For extremely thin Salmon layers (h0=0.1)
  Nt = ceil(tmax/dt)+1;       
  
  %%% Set viscosityy  
  A2 = 0;
  A4 = 0;
  A4smag = 4;

  %%% Salmon layer thickness  
  %%% Chosen to satisfy Salmon's (2002) stability criterion, approximately.
  switch (config)
    case 'wind'
      if (Ny==256)
        h0 = 1.5;
      else
        h0 = 3;
      end
%       h0 = 1;
%       h0 = 0.5;
%       h0 = 0.1;
    case 'rand'
      h0 = 5 * 128/Ny;
  end
%   h0 = 3;  %%% For many-layer config
  hsml = 50;
  hbbl = 50;
  hmin_surf = 12; %%% Should be >> h0 to avoid forcing of very thin layers
%   hmin_surf = 6;  %%% For many-layer config
  hmin_bot = 0;
        
         
  %%% Define parameters 
  PARAMS = addParameter(PARAMS,'Nlay',Nlay,PARM_INT);
  PARAMS = addParameter(PARAMS,'Nx',Nx,PARM_INT);
  PARAMS = addParameter(PARAMS,'Ny',Ny,PARM_INT);
  PARAMS = addParameter(PARAMS,'dt',dt,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefrequency',savefreq,PARM_REALF);  
  PARAMS = addParameter(PARAMS,'savefreqEZ',savefreqEZ,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreqAvg',savefreqAvg,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreqUMom',savefreqUMom,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreqVMom',savefreqVMom,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreqThic',savefreqThic,PARM_REALF);  
  PARAMS = addParameter(PARAMS,'tmax',tmax,PARM_REALF);
  PARAMS = addParameter(PARAMS,'restart',restart,PARM_INT);
  PARAMS = addParameter(PARAMS,'startIdx',startIdx,PARM_INT);
  PARAMS = addParameter(PARAMS,'Ly',Ly,PARM_REALF);
  PARAMS = addParameter(PARAMS,'h0',h0,PARM_REALF);
  PARAMS = addParameter(PARAMS,'hsml',hsml,PARM_REALF);
  PARAMS = addParameter(PARAMS,'hbbl',hbbl,PARM_REALF);
  PARAMS = addParameter(PARAMS,'hmin_surf',hmin_surf,PARM_REALF);
  PARAMS = addParameter(PARAMS,'hmin_bot',hmin_bot,PARM_REALF);
  PARAMS = addParameter(PARAMS,'A2',A2,PARM_REALF);
  PARAMS = addParameter(PARAMS,'A4',A4,PARM_REALF);
  PARAMS = addParameter(PARAMS,'A4smag',A4smag,PARM_REALF);
  PARAMS = addParameter(PARAMS,'linDragCoeff',rb,PARM_REALF);
  PARAMS = addParameter(PARAMS,'quadDragCoeff',Cd,PARM_REALF);
  PARAMS = addParameter(PARAMS,'quadDragSurf',Cd_surf,PARM_REALF);
  PARAMS = addParameter(PARAMS,'use_MG',use_MG,PARM_INT);
  PARAMS = addParameter(PARAMS,'useRL',useRL,PARM_INT);
  PARAMS = addParameter(PARAMS,'useWallEW',0,PARM_INT);
  PARAMS = addParameter(PARAMS,'useWallNS',1,PARM_INT);
  PARAMS = addParameter(PARAMS,'tol',tol,PARM_REALE);
  PARAMS = addParameter(PARAMS,'SOR_rp_max',SOR_rp_max,PARM_REALF);
  PARAMS = addParameter(PARAMS,'SOR_rp_min',SOR_rp_min,PARM_REALF);   
  PARAMS = addParameter(PARAMS,'SOR_opt_freq',SOR_opt_freq,PARM_INT);   
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% BACKGROUND ROTATION %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  
  %%% Background rotation rate  
  Omegaz = 0.5 * (f0 + beta_x*(XX_q-Lx/2) + beta_y*(YY_q-Ly/2));
  
  fignum = fignum+1;
  figure(fignum);
  pcolor(XX_q,YY_q,2*Omegaz);
  shading interp;
  colorbar;
  

  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WIND STRESS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Coast-localized wind
  Lwind = 200*m1km;
  taux = tau0*sin(pi*YY_u/Lwind).^2 / rho0;  
  taux(YY_u>Lwind) = 0;

  %%% Realistic wind shape with coastal drop-off
%   Lwind = 75*m1km;
%   taux = tau0*tanh(YY_u/Lwind).^2 / rho0;
  
  if (go_west)
    taux = -taux;
  end
  
  fignum = fignum+1;
  figure(fignum);
  plot(yy_h,taux(1,:));
  title('Wind stress');
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% ALONG-SHORE PRESSURE GRADIENT %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Zonal baroropic forcing  
  Fbaro_tot = Fbaro0*sin(pi*YY_u/Lwind).^2 / rho0;
  Fbaro_tot(YY_u>Lwind) = 0;  
  Fbaro = Fbaro_tot ./ repmat(mean(-etab,1),[Nx 1]);
  
  fignum = fignum+1;
  figure(fignum);
  plot(yy_h,Fbaro(1,:));
  title('Barotropic forcing, x = 0');
  xlabel('y');
  ylabel('m/s^2');

  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RELAXATION SETUP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Initialize with zero flow uniform layer thicknesses
  uRelax = zeros(Nlay,Nx,Ny);
  vRelax = zeros(Nlay,Nx,Ny);  
  eRelax = repmat(reshape(E0(1:Nlay),[Nlay 1 1]),[1 Nx Ny]);
  hRelax = repmat(reshape(H0,[Nlay 1 1]),[1 Nx Ny]);
  
  %%% Relaxation time scale
  relTime = -ones(Nx,Ny);    
  relTime(YY_h>Ly-Lrelax) = 1 ./ (1-(Ly-YY_h(YY_h>Ly-Lrelax))/Lrelax);
  eTime = tRelax_h * repmat(reshape(relTime,[1 Nx Ny]),[Nlay 1 1]);
  uTime = tRelax_u * repmat(reshape(relTime,[1 Nx Ny]),[Nlay 1 1]);
  vTime = tRelax_v * repmat(reshape(relTime,[1 Nx Ny]),[Nlay 1 1]);
  hTime = tRelax_h * relTime;
  
  fignum = fignum+1;
  figure(fignum);
  plot(YY_h(1,:),squeeze(eRelax(:,Nx/2,:)));
  title('Relaxation layer elevations');
  xlabel('y')
  xlabel('z');
  
  fignum = fignum+1; 
  figure(fignum);
  plot(YY_h(1,:),squeeze(eTime(:,Nx/2,:)));
  title('Relaxation time scales');
  xlabel('y');
  ylabel('(s)')

    
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% INITIAL CONDITIONS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  %%% Initialize with zero flow
  uu_init = uRelax*0;
  vv_init = vRelax*0; 
  
  %%% Initialize Montgomery potential with reference values  
  MM_init = zeros(Nlay,Nx,Ny);
  for k=2:Nlay
    MM_init(k,:,:) = MM_init(k-1,:,:) + gg(k)*E0(k);
  end 
     
  %%% Compute and regularize interface depths prior to calculation of
  %%% Salmon depths
  hh_init = repmat(reshape(H0,[Nlay 1 1]),[1 Nx Ny]);
  ee_init = zeros(Nlay,Nx,Ny);
  ee_init(2:Nlay+1,:,:) = -cumsum(hh_init);
  ee_init(Nlay+1,:,:) = etab;
  for i=1:Nx
    for j=1:Ny
      for k=Nlay:-1:2
        if (ee_init(k,i,j)<ee_init(k+1,i,j))
          ee_init(k,i,j) = ee_init(k+1,i,j) + h0;
        end
      end
    end
  end
  hh_init = ee_init(1:Nlay,:,:) - ee_init(2:Nlay+1,:,:);
  
  %%% Adjust layer thicknesses to account for incrops
  hh_init = solveSalmonThicknesses (hh_init,etab,zeros(Nx,Ny),MM_init,gg,h0);   
  
  %%% Plot thicknesses
  fignum = fignum+1;
  figure(fignum);
  plot(xx_h,hh_init(:,:,1));
  title('Initial layer thicknesses');
  
  %%% Plot layer elevations
  ee_init = zeros(Nlay+1,Nx,Ny);
  ee_init(2:Nlay+1,:,:) = -cumsum(hh_init);  
  
  fignum = fignum+1;
  figure(fignum);
  plot(xx_h,ee_init(:,:,1));
  title('Initial layer interfaces');
  xlabel('x');
  ylabel('z');
  
  fignum = fignum+1;
  figure(fignum);
  plot(yy_h,squeeze(ee_init(:,Nx/2,:)));
  title('Initial layer interfaces');
  xlabel('y');
  ylabel('z');
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RANDOM FORCING SETUP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %%% Generate forcing mask in spectral space
  k = [0:1:Nx/2-1,-Nx/2:1:-1]; %%% Zonal wavenumber
  K_xk = 2*pi.*(k)./Lx;
  l = [0:1:Ny/2-1,-Ny/2:1:-1]; %%% Meridional wavenumber
  K_yl = 2*pi.*(l)./Ly;
  [K_ykl,K_xkl] = meshgrid(K_yl, K_xk);
  K = sqrt(K_xkl.^2 + K_ykl.^2); %%% Absolute wavenumber 
  K_0 = 2*pi/lambdaK; %%% Most strongly forced wavenumber
  W = K_0/8; %%% Exponential width of energy band in wavenumber space   
  RF_mask_fft = exp(-((K-K_0)/W).^2);
  RF_mask_fft = repmat(reshape(RF_mask_fft,[1 Nx Ny]),[Nlay 1 1]);
    
  %%% Generate forcing mask in real space
  idx_south_zero = (yy_q(1:Ny) <= Ymin_RF);
  idx_south_cos = (yy_q(1:Ny) <= Ymin_RF+Lramp_RF) & (yy_q(1:Ny) >= Ymin_RF);
  idx_north_cos = (yy_q(1:Ny) >= Ymax_RF-Lramp_RF) & (yy_q(1:Ny) <= Ymax_RF);
  idx_north_zero = (yy_q(1:Ny) >= Ymax_RF);
  RF_mask_q = ones(Nlay,Nx,Ny);
  if (~isempty(idx_south_zero))
    RF_mask_q(:,:,idx_south_zero) = 0;
  end
  if (~isempty(idx_south_cos))    
    RF_mask_q(:,:,idx_south_cos) = RF_mask_q(:,:,idx_south_cos) .* repmat(reshape(0.5.*(1+cos(pi*(YY_q(1:Nx,idx_south_cos)-Ymin_RF-Lramp_RF)/Lramp_RF)),[1 Nx sum(idx_south_cos)]),[Nlay 1 1]);
  end
  if (~isempty(idx_north_cos))
    RF_mask_q(:,:,idx_north_cos) = RF_mask_q(:,:,idx_north_cos) .* repmat(reshape(0.5.*(1+cos(pi*(YY_q(1:Nx,idx_north_cos)-Ymax_RF-Lramp_RF)/Lramp_RF)),[1 Nx sum(idx_north_cos)]),[Nlay 1 1]);
  end
  if (~isempty(idx_north_zero))
    RF_mask_q(:,:,idx_north_zero) = 0;
  end 
  
  %%% Modify parameters depending on whether we want barotropic or
  %%% baroclinic forcing
  if (Nlay == 1)
    useThicWeightedRF =  true; %%% Apply random forcing to momentum equation
  else
    %%% Apply the forcing to the velocity equation and divide the forcing 
    %%% mask in real space by the water column thickness at each point.
    if (useBaroclinicForcing)
      %%% Apply forcing to near-surface layers only
      useThicWeightedRF = false;
      RF_mask_q(Nlay,:,:) = 0;      
      RF_mask_q = RF_mask_q / H;
    else    
      %%% Apply forcing to full water column depth
      %%% This ensures that the depth-integrated forcing, effectively
      %%% sum_k h_ijk * mask_ij * F0 * curl Psi_ij, has a magnitude of F0.
      useThicWeightedRF = false;      
      for i=1:Nx
        for j=1:Ny
          RF_mask_q(:,i,j) = RF_mask_q(:,i,j) / H * (-H/etab(i,j)).^0;
        end
      end
    end
  end
  
  
  %%% Make a plot
  fignum = fignum+1;
  figure(fignum);
  plot(yy_q(1:Ny),squeeze(RF_mask_q(:,1,:)));
  xlabel('y');
  title('Random forcing factor, x = 0');


  %%% Add scalar parameters
  PARAMS = addParameter(PARAMS,'useRandomForcing',useRandomForcing,PARM_INT);
  PARAMS = addParameter(PARAMS,'useThicWeightedRF',useThicWeightedRF,PARM_INT);
  PARAMS = addParameter(PARAMS,'RF_F0_rot',RF_F0_rot,PARM_REALF);
  PARAMS = addParameter(PARAMS,'RF_tau',RF_tau,PARM_REALF);

  %%% Add matrix parameters
  RF_fftMaskFile = 'RF_fftMaskFile.dat';
  writeDataFile(fullfile(local_run_dir,RF_fftMaskFile),RF_mask_fft);
  PARAMS = addParameter(PARAMS,'RF_fftMaskFile',RF_fftMaskFile,PARM_STR);  
  RF_rotMaskFile = 'RF_rotMaskFile.dat';
  writeDataFile(fullfile(local_run_dir,RF_rotMaskFile),RF_mask_q);
  PARAMS = addParameter(PARAMS,'RF_rotMaskFile',RF_rotMaskFile,PARM_STR);
  
  
  
  
  
  
  
  

  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CREATE INPUT FILES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  
  %%% Initial h 
  hInitFile = 'hInit.dat';
  writeDataFile(fullfile(local_run_dir,hInitFile),hh_init);
  PARAMS = addParameter(PARAMS,'hInitFile',hInitFile,PARM_STR);  
  
  %%% Initial u  
  uInitFile = 'uInit.dat';
  writeDataFile(fullfile(local_run_dir,uInitFile),uu_init);
  PARAMS = addParameter(PARAMS,'uInitFile',uInitFile,PARM_STR); 
  
  %%% Initial v
  vInitFile = 'vInit.dat';
  writeDataFile(fullfile(local_run_dir,vInitFile),vv_init);
  PARAMS = addParameter(PARAMS,'vInitFile',vInitFile,PARM_STR);  
  
  %%% Bathymetry
  etabFile = 'etab.dat';          
  writeDataFile(fullfile(local_run_dir,etabFile),etab);
  PARAMS = addParameter(PARAMS,'hbFile',etabFile,PARM_STR);
  
  %%% Background rotation
  OmegazFile = 'Omegaz.dat';          
  writeDataFile(fullfile(local_run_dir,OmegazFile),Omegaz);
  PARAMS = addParameter(PARAMS,'OmegazFile',OmegazFile,PARM_STR);
  
  %%% Reduced gravity
  gFile = 'gg.dat';          
  writeDataFile(fullfile(local_run_dir,gFile),gg);
  PARAMS = addParameter(PARAMS,'gFile',gFile,PARM_STR);

  %%% Wind stress
  tauxFile = 'taux.dat';          
  writeDataFile(fullfile(local_run_dir,tauxFile),taux);
  PARAMS = addParameter(PARAMS,'tauxFile',tauxFile,PARM_STR);
  
  %%% Along-shore pressure gradient
  FbaroXfile = 'Fbaro_x.dat';          
  writeDataFile(fullfile(local_run_dir,FbaroXfile),Fbaro);
  PARAMS = addParameter(PARAMS,'FbaroXFile',FbaroXfile,PARM_STR);
  
  %%% Relaxation values for u
  uRelaxFile = 'uRelax.dat';
  writeDataFile(fullfile(local_run_dir,uRelaxFile),uRelax);
  PARAMS = addParameter(PARAMS,'uRelaxFile',uRelaxFile,PARM_STR);  
  
  %%% Relaxation timescale for u  
  uTimeFile = 'uTime.dat';
  writeDataFile(fullfile(local_run_dir,uTimeFile),uTime);
  PARAMS = addParameter(PARAMS,'uTimeFile',uTimeFile,PARM_STR); 
  
  %%% Relaxation values for v
  vRelaxFile = 'vRelax.dat';
  writeDataFile(fullfile(local_run_dir,vRelaxFile),vRelax);
  PARAMS = addParameter(PARAMS,'vRelaxFile',vRelaxFile,PARM_STR);  
  
  %%% Relaxation timescale for v  
  vTimeFile = 'vTime.dat';
  writeDataFile(fullfile(local_run_dir,vTimeFile),vTime);
  PARAMS = addParameter(PARAMS,'vTimeFile',vTimeFile,PARM_STR); 
  
  %%% Relaxation values for eta
  eRelaxFile = 'eRelax.dat';
  writeDataFile(fullfile(local_run_dir,eRelaxFile),eRelax);
  PARAMS = addParameter(PARAMS,'eRelaxFile',eRelaxFile,PARM_STR);  
  
  %%% Relaxation timescale for eta
  eTimeFile = 'eTime.dat';
  writeDataFile(fullfile(local_run_dir,eTimeFile),eTime);
  PARAMS = addParameter(PARAMS,'eTimeFile',eTimeFile,PARM_STR);  
  
%   %%% Relaxation values for h
%   hRelaxFile = 'hRelax.dat';
%   writeDataFile(fullfile(local_run_dir,hRelaxFile),hRelax);
%   PARAMS = addParameter(PARAMS,'hRelaxFile',hRelaxFile,PARM_STR);  
%   
%   %%% Relaxation timescale for h
%   hTimeFile = 'hTime.dat';
%   writeDataFile(fullfile(local_run_dir,hTimeFile),hTime);
%   PARAMS = addParameter(PARAMS,'hTimeFile',hTimeFile,PARM_STR);  

  %%% Create the input parameter file
  writeParamFile(pfname,PARAMS);      

  %%% Create a run script
  createRunScript (  local_home_dir, ...
                     run_name, ...
                     model_code_dir, ...
                     exec_name, ...
                     use_intel, ...
                     use_pbs, ... 
                     use_cluster, ...
                     uname, ...
                     cluster_addr, ...
                     cluster_home_dir);

end