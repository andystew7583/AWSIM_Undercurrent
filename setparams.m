%%%
%%% setparams.m
%%%
%%% Configures a shelf/slope channel with wind and buoyancy forcing.
%%% 
function setparams (local_home_dir,run_name)

  %%% Toggles number of layers
  use_three_layer = false; 

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
  
  %%% Physical parameters
  rho0 = 1000;                  %%% Reference density 
  f0 = -1e-4;                   %%% Coriolis parameter
  beta = 1e-11;               %%% Coriolis parameter gradient  
  Ly = 250*m1km;                %%% Domain length.   
  Lx = 2*Ly;
  if (use_three_layer)
    geff = [g 0.5e-2 0.5e-2];   %%% Reduced gravity at layer interface 
  else
    geff = [g 0.5e-2];          %%% Reduced gravity at layer interface 
  end     
  H = 4000;                     %%% Ocean depth  
  if (use_three_layer)    
    H1 = 200;                   %%% Initial upper layer thickness
    H2 = 200;                   %%% Initial middle layer thickness
    H0 = [H1 H2 H-H1-H2];       %%% Initial layer thicknesses - used for wave speed calculation  
  else
    H1 = 250;                   %%% Initial upper layer thickness    
    H0 = [H1 H-H1];             %%% Initial layer thicknesses - used for wave speed calculation  
  end
  E0 = 0.01;                    %%% Initial EKE density
  rb = 0;                    %%% Linear bottom drag
  Cd = 2e-3;                %%% Quadratic bottom drag
  Cd_surf = 0e-3;           %%% Surface drag coefficient
  tau0 = 0;                  %%% Wind stress maximum    
  go_west = true;              %%% Set true for westward flow - reverses wind stress and initial flow  
  if (use_three_layer)
    eta_north = [-200 -400];     %%% Relaxation layer depths at northern boundary
  else
    eta_north = -200;
    eta_south = -350;
    eta_south_dense = -50;
  end  
  deta1 = (eta_north-eta_south)*geff(2)/geff(1);                    %%% Initial SSH change across the channel (or equivalent surface pressure change)  
  if (go_west)                  %%% Relaxation time scales
    tRelax_u = -1;
    tRelax_v = -1;
    tRelax_h = 7*t1day;                
  else
    tRelax_u = -1;
    tRelax_v = -1;
    tRelax_h = 7*t1day;         
  end
  Lrelax = 25*m1km;            %%% Meridional width of nudging zone
  Xrelax = Lx/2;               %%% Location of buoyancy forcing on shelf
  Xrelax_dense = 0;           
  Wrelax = 25*m1km;               %%% Zonal width of shelf forcing
  
  %%% Temporal parameters  
  tmax = 10*t1year;
%   tmax = 100*t1day;
  savefreq = 5*t1day; 
%   savefreq = 400; 
  savefreqEZ = 1*t1day;
%   savefreqEZ = 0.1*t1day;
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
  if (use_three_layer)
    Nlay = 3;
  else
    Nlay = 2;
  end
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
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% BOTTOM TOPOGRAPHY %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %%% Put canyons in continental shelf
  H_shelf = 500; %%% Continental shelf water column thickness  
  W_canyons = 60*m1km;
  H_canyons = 300;
  X_canyons = [Lx/4 3*Lx/4];  
  Z_shelf = -H_shelf * ones(Nx,Ny); %%% Continental shelf elevation
  for m=1:length(X_canyons)
    idx = (XX_h > X_canyons(m)-W_canyons/2) & (XX_h < X_canyons(m)+W_canyons/2);
    Z_shelf(idx) = Z_shelf(idx) - H_canyons*0.5*(1+cos(2*pi*(XX_h(idx)-X_canyons(m))/(W_canyons))); 
  end
  
  %%% Create continental slope
%   Y_slope = 100*m1km; %%% Latitude of center of continental slope
%   W_slope = 20*m1km;
%   gam_h = 0.004;
%   Y_scaled = (YY_h-Y_slope)/W_slope;     
%   Z_slope = 0.5*(Z_shelf-H);
%   H_slope = Z_shelf + H;
%   etab = Z_slope + H_slope .* (0.25*sqrt((1-Y_scaled).^2 + 4*gam_h*Y_scaled.^2)-0.25*sqrt((1+Y_scaled).^2 +  4*gam_h*Y_scaled.^2)) / (1+4*gam_h)^(-1/2);
  
  %%% Exponentials and linear slope with matched gradients
  smax = 1.5e-1; %%% Max slope steepness
  Ws_shelf = 2*m1km; %%% Reference slope exponential width for continental shelf
  Ws_deep = 12.5*m1km; %%% Reference slope exponential width for deep ocean
  Ys_shelf = 80*m1km; %%% Latitude of center of continental slope  
  
  %%% Matrices defining spatially-dependent geometric parameters
  Z_bot = -H * ones(Nx,Ny); %%% Flat reference ocean floor depth    
  Ws_upper = Ws_shelf * ones(Nx,Ny); %%% Exponential width of upper slope
  Ws_lower = Ws_deep * ones(Nx,Ny); %%% Exponential width of lower slope, reverts to Ws_shelf where ridge elevation is high
  Delta_upper = Ws_upper * smax; %%% Exponential elevation change of upper slope
  Delta_lower = Ws_lower * smax; %%% Exponential elevation change of lower slope
  Zs_upper = (Z_shelf - Delta_upper); %%% Elevation of upper/mid slope join
  Zs_lower = (Z_bot + Delta_lower); %%% Elevation of lower/mid slope join  
  Ys_upper = Ys_shelf; %%% Latitude of upper/mid slope join
  Ys_lower = Ys_upper + ((Zs_upper-Zs_lower)/smax); %%% Latitude of lower/mid slope join
  
  %%% Construct the bathymetry by joining exponential upper and lower
  %%% slopes with a linear mid-slope
  etab = zeros(Nx,Ny);
  idx_upper = YY_h < Ys_upper;
  idx_mid = (YY_h >= Ys_upper) & (YY_h < Ys_lower);
  idx_lower = YY_h > Ys_lower;    
  etab_upper = Z_shelf - Delta_upper .* exp((YY_h-Ys_upper)./Ws_upper);
  etab_mid = Zs_upper - smax.*(YY_h-Ys_upper);
  etab_lower = Z_bot + Delta_lower .* exp(-(YY_h-Ys_lower)./Ws_lower);
  etab(idx_upper) = etab_upper(idx_upper);
  etab(idx_mid) = etab_mid(idx_mid);
  etab(idx_lower) = etab_lower(idx_lower);
       
  %%% Plot topography
  figure(10);
  surf(XX_h,YY_h,etab);
  shading interp;
  colorbar;
%   plot(yy_h,etab(1,:))

  figure(11)
  plot(yy_h,etab(Nx/2,:));
  hold on;
  plot(yy_h,etab(Nx/4,:));
  hold off;
  
  figure(12)
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
  c = calcWaveSpeed(H0,geff,useRL)
  Umax = 3;  
  dt = 0.25*d/(c+Umax)  
  Nt = ceil(tmax/dt)+1;      
  H0
  geff
  
  %%% Set viscosityy  
%   tau_diff = 2*t1day;
%   tau_diff = 10*t1day/sqrt(alpha);  
%   A2 = d^2/tau_diff;
  A2 = 0;
%   A4 = d^4/tau_diff
%   A4 = 0.01*d^3*Umax
%   A4 = 0;
%   A4grid = 0.1; %%% Grid viscosity (must be < 1, typically much less)
%   A4 = 0.25*0.125*d^4/dt * A4grid
  A4 = 0;
  A4smag = 4;

 %%% Salmon layer thickness  
  %%% Chosen to satisfy Salmon's (2002) stability criterion, approximately.
  %%% The prefactor has been chosen empirically.
  h0 = 5;
  hsml = 30;
  hbbl = 30;
  hmin_surf = 10;
  hmin_bot = 0;
        
         
  %%% Define parameters 
  PARAMS = addParameter(PARAMS,'Nlay',Nlay,PARM_INT);
  PARAMS = addParameter(PARAMS,'Nx',Nx,PARM_INT);
  PARAMS = addParameter(PARAMS,'Ny',Ny,PARM_INT);
  PARAMS = addParameter(PARAMS,'dt',dt,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefrequency',savefreq,PARM_REALF);
  PARAMS = addParameter(PARAMS,'savefreqEZ',savefreqEZ,PARM_REALF);
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
  Omegaz = 0.5* (f0*ones(Nx+1,Ny+1) + beta*(YY_q-(Ly/2)));
  


  

  %%%%%%%%%%%%%%%%%%%%%%%
  %%%%% WIND STRESS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Zonal wind stress
  Lwind = 200*m1km;
  taux = tau0*sin(pi*YY_u/Lwind).^2 / rho0;
  taux(YY_u>Lwind) = 0;
  if (go_west)
    taux = -taux;
  end
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% INITIAL CONDITIONS %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Set sea surface height
  Rd = c/abs(f0)
  lambdaK = 4*Rd;  
  eta1 = (f0/g)*genRandIC(lambdaK,E0,Nx,Ny,Ly);
  eta1 = (1 - exp(-(YY_h./(lambdaK)).^2)) .* (1 - exp(-((YY_h-Ly)./(lambdaK)).^2)) .* eta1;  
  eta1 = eta1 + deta1*(YY_h-Ly/2)/Ly;
  if (go_west)
    eta1 = -eta1;
  end
  if (use_three_layer)
    eta2 = -H1 - (geff(1)/geff(2))*eta1;
    eta3 = (-H1-H2)*ones(Nx,Ny);
    h1 = -eta2;
    h2 = eta2-eta3;
    h3 = eta3-etab;
  else
    eta2 = -H1 - (geff(1)/geff(2))*eta1;    
    h1 = -eta2;
    h2 = eta2-etab;    
  end
 
  %%% Plot layer interfaces
  figure(1);
  contourf(XX_h,YY_h,eta1);
  colorbar;  
  figure(2);
  contourf(XX_h,YY_h,eta2);
  colorbar;    
  
  %%% Set geostrophic along-slope velocity
  eta1_mat1 = circshift(eta1, [0,-1]);  
  eta1_mat2 = circshift(eta1, [1,-1]);
  eta1_mat3 = circshift(eta1, [1,0]);
  eta1_mat4 = circshift(eta1, [1,1]);
  eta1_mat5 = circshift(eta1, [0,1]); 
  u1 = -(g/f0*(1/d)) .* ((1/4).*(eta1_mat1+eta1_mat2+eta1_mat3+eta1) - (1/4).*(eta1+eta1_mat3+eta1_mat4+eta1_mat5));
  u2 = 0*u1;
  u3 = 0*u1;
  
  %%% Quick fix for near-boundary points
  u1(:,1) = u1(:,2);    
  u1(:,Ny) = u1(:,Ny-1);
  u2(:,1) = u2(:,2);
  u2(:,Ny) = u2(:,Ny-1);
  u3(:,1) = u3(:,2);
  u3(:,Ny) = u3(:,Ny-1);
  
  %%% Plot along-slope velocity
  figure(3);
  contourf(XX_u,YY_u,u1);
  colorbar;
  
  %%% Set geostrophic cross-slope velocity
  eta1_mat6 = circshift(eta1, [-1,1]);
  eta1_mat7 = circshift(eta1, [-1,0]);  
  v1 = -(g/f0*(1/d)) .* ((1/4).*(eta1_mat3+eta1_mat4+eta1_mat5+eta1) - (1/4).*(eta1+eta1_mat5+eta1_mat6+eta1_mat7));
  v2 = 0*v1;
  v3 = 0*v1;
  
  %%% Plot cross-slope velocity
  figure(4);
  contourf(XX_v,YY_v,v1);     
  colorbar;
         
  %%% Calculate boundary vorticities     
  u1_mat5 = circshift(u1, [0,1]);
  v1_mat3 = circshift(v1, [1,0]);
  zeta1 = zeros(Nx+1,Ny+1);
  zeta1(1:Nx,1:Ny) = (v1-v1_mat3)./d - (u1-u1_mat5)./d;  
  zeta1(:,1) = 0;
  zeta1(:,end) = 0;
  
  %%% Plot vorticity
  figure(5);  
  contourf(XX_q/1000,YY_q/1000,zeta1/abs(f0));   
  colorbar;
  colormap redblue;
  
  %%% Calculate kinetic energy of initial state
  u_mat7 = circshift(u1, [-1,0]);
  v_mat1 = circshift(v1, [0,-1]);
  %KE = 1/2*d*d.*(zeros(Nx, Ny)-etab).*((u.^2+u_mat7.^2)./2 + (v.^2+v_mat1.^2)./2);
  KE = (1/2).*((u1.^2+u_mat7.^2)./2 + (v1.^2+v_mat1.^2)./2);
  meanE = (1/(Nx*Ny))*sum(sum(KE))
  
  %%% Plot fourier transform of Kinetic Energy
  vfft = fft2(v1);
  ufft = fft2(u1);   
  KEfft = vfft.^2 + ufft.^2;  
  figure(6);
  [p,q] =contourf(real(KEfft),20);
  colorbar;     
  
  %%% Apply the same initial velocity profiles in each layer
  uu = zeros(Nlay,Nx,Ny);
  uu(1,:,:) = u1;
  uu(2,:,:) = u2;
  if (use_three_layer)
    uu(3,:,:) = u3;
  end
  vv = zeros(Nlay,Nx,Ny);
  vv(1,:,:) = v1;
  vv(2,:,:) = v2;
  if (use_three_layer)
    vv(3,:,:) = v3;
  end
  hh = zeros(Nlay,Nx,Ny);
  hh(1,:,:) = h1;
  hh(2,:,:) = h2;  
  if (use_three_layer)
    hh(3,:,:) = h3;  
  end
   
  %%% Plot initial layer thickness
  figure(7);
  plot(xx_h,squeeze(hh(2,:,round(Ny/2))));

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% RELAXATION SETUP %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% Relaxation targets for layer thicknesses
  h1Relax = 0*h1;
  h2Relax = 0*h2;   
  h1Relax(YY_h<Ly/2) = - eta_south(1);
  h1Relax(YY_h>=Ly/2) = - eta_north(1);  
  if (Nlay == 2)
    h2Relax(YY_h<Ly/2) = eta_south(1) - etab(YY_h<Ly/2);
    h2Relax(YY_h>=Ly/2) = eta_north(1) - etab(YY_h>=Ly/2);      
    h2Relax((YY_h<Ly/2) & ((XX_h>3*Lx/4) | (XX_h<Lx/4))) = eta_south_dense(1) - etab((YY_h<Ly/2) & ((XX_h>3*Lx/4) | (XX_h<Lx/4)));
  end
  if (Nlay == 3)
    h3Relax = 0*h3; 
    h2Relax(YY_h<Ly/2) = eta_south(1) - eta_south(2);
    h2Relax(YY_h>=Ly/2) = eta_north(1) - eta_north(2);  
    h3Relax(YY_h<Ly/2) = eta_south(2) - etab(YY_h<Ly/2);
    h3Relax(YY_h>=Ly/2) = eta_north(2) - etab(YY_h>=Ly/2);
  end  
  
  %%% Relaxation time scale
%   hTime = -ones(Nx,Ny);  
%   hTime(YY_h<Lrelax) = tRelax_h ./ (1-YY_h(YY_h<Lrelax)/Lrelax);
%   hTime = tRelax_h ./ ( exp(-(YY_h/(Lrelax/2)).^2) .* exp(-((XX_h-Xrelax)/Wrelax).^2) );
  XX_dense = XX_h;
  XX_dense(XX_dense>Lx/2) = XX_dense(XX_dense>Lx/2) - Lx;
  hTime = tRelax_h ./ ( exp(-(YY_h/(Lrelax/2)).^2) .* (exp(-((XX_h-Xrelax)/Wrelax).^2) + exp(-((XX_dense-Xrelax_dense)/Wrelax).^2)) );
  hTime(YY_h>Lrelax) = -1;
  hTime(YY_h>Ly-Lrelax) = tRelax_h ./ (1-(Ly-YY_h(YY_h>Ly-Lrelax))/Lrelax);
  
  figure(99)
  plot(YY_h(1,:),h1Relax(Nx/2,:));
  hold on;
  plot(YY_h(1,:),h2Relax(Nx/2,:));
  hold off;
  
  figure(98)
  plot(XX_h(:,1),h2Relax(:,1));
  
  figure(100);
  semilogy(YY_h(Nx/2,:),hTime(Nx/2,:)/t1year);
  
  figure(101);
  semilogy(XX_h(:,1),hTime(:,1)/t1year);
  
  %%% Create input matrices
  hRelax = zeros(Nlay,Nx,Ny);
  hRelax(1,:,:) = h1Relax;
  hRelax(2,:,:) = h2Relax;
  if (Nlay == 3)    
    hRelax(3,:,:) = h3Relax;
  end
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%% CREATE INPUT FILES %%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
  
  %%% Initial h 
  hInitFile = 'hInit.dat';
  writeDataFile(fullfile(local_run_dir,hInitFile),hh);
  PARAMS = addParameter(PARAMS,'hInitFile',hInitFile,PARM_STR);  
  
  %%% Initial u  
  uInitFile = 'uInit.dat';
  writeDataFile(fullfile(local_run_dir,uInitFile),uu);
  PARAMS = addParameter(PARAMS,'uInitFile',uInitFile,PARM_STR); 
  %%% Initial v
  vInitFile = 'vInit.dat';
  writeDataFile(fullfile(local_run_dir,vInitFile),vv);
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
  writeDataFile(fullfile(local_run_dir,gFile),geff);
  PARAMS = addParameter(PARAMS,'gFile',gFile,PARM_STR);

  %%% Wind stress
  tauxFile = 'taux.dat';          
  writeDataFile(fullfile(local_run_dir,tauxFile),taux);
  PARAMS = addParameter(PARAMS,'tauxFile',tauxFile,PARM_STR);
  
  %%% Relaxation values for h 
  hRelaxFile = 'hRelax.dat';
  writeDataFile(fullfile(local_run_dir,hRelaxFile),hRelax);
  PARAMS = addParameter(PARAMS,'hRelaxFile',hRelaxFile,PARM_STR);  
  
  %%% Relaxation timescale for h  
  hTimeFile = 'hTime.dat';
  writeDataFile(fullfile(local_run_dir,hTimeFile),hTime);
  PARAMS = addParameter(PARAMS,'hTimeFile',hTimeFile,PARM_STR);  

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