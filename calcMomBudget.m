%%%
%%% calcMomBudget.m
%%%
%%% Calculates diagnostics from time-averaged momentum budget for one of our
%%% undercurrent runs.
%%%
function calcMomBudget (local_home_dir,run_name,tmax,uc_layidx)

  %%% Define coordinate system for analysis
  coord_name = 'bathy';
  % coord_name = 'psi';

  %%% Load experiment parameters
  dirpath = fullfile(local_home_dir,run_name);
  loadParams;
  rho0 = 1000;

  %%% Define coordinate grid
  switch (coord_name)
    case 'bathy'
      dd = [110:10:3500]';   
    case 'psi'
      dd = [0.01:0.01:0.5 0.55:0.05:5];
  end
  Nd = length(dd); 
  
  %%% Define averaging period
%   tmax = tmax + 0.05*t1year;
  tmax = tmax - 2*t1year+ 0.05*t1year;
%   tmin = tmax - 10*t1year;
  tmin = tmax - 8*t1year;
   
  %%% Average u-momentum diagnostics  
%   avg_iter_start = n0_avg_hu;
%   avg_num_iters = N_avg_hu;
%   avg_start_time = startTime;
  avg_iter_start = 0;                    
  avg_num_iters = N_avg_hu+n0_avg_hu;
  avg_start_time = 0;
  UMom_PVadvection = rho0*do_avg(dirpath,OUTN_UMOM_Q,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hu,tmin,tmax,avg_start_time);
  UMom_Montgomery = rho0*do_avg(dirpath,OUTN_UMOM_GRADM,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hu,tmin,tmax,avg_start_time);
  UMom_KEgradient = rho0*do_avg(dirpath,OUTN_UMOM_GRADKE,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hu,tmin,tmax,avg_start_time);
  UMom_dhdt = rho0*do_avg(dirpath,OUTN_UMOM_DHDT,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hu,tmin,tmax,avg_start_time);
  UMom_hypervisc = rho0*do_avg(dirpath,OUTN_UMOM_A4,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hu,tmin,tmax,avg_start_time);
  UMom_quadBotDrag = rho0*do_avg(dirpath,OUTN_UMOM_CDBOT,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hu,tmin,tmax,avg_start_time);
  UMom_windStress = rho0*do_avg(dirpath,OUTN_UMOM_WIND,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hu,tmin,tmax,avg_start_time);
  UMom_relaxation = rho0*do_avg(dirpath,OUTN_UMOM_RELAX,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hu,tmin,tmax,avg_start_time);
  UMom_baroForcing = rho0*do_avg(dirpath,OUTN_UMOM_FBARO,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hu,tmin,tmax,avg_start_time);
  UMom_randomForcing = rho0*do_avg(dirpath,OUTN_UMOM_RAND,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hu,tmin,tmax,avg_start_time);
  UMom_diapycnal = rho0*do_avg(dirpath,OUTN_UMOM_WDIA,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hu,tmin,tmax,avg_start_time); 

  %%% Average v-momentum diagnostics  
%   avg_iter_start = n0_avg_hv;
%   avg_num_iters = N_avg_hv;
%   avg_start_time = startTime;
  avg_iter_start = 0;                    
  avg_num_iters = N_avg_hu+n0_avg_hu;
  avg_start_time = 0;
  VMom_PVadvection = rho0*do_avg(dirpath,OUTN_VMOM_Q,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hv,tmin,tmax,avg_start_time);
  VMom_Montgomery = rho0*do_avg(dirpath,OUTN_VMOM_GRADM,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hv,tmin,tmax,avg_start_time);
  VMom_KEgradient = rho0*do_avg(dirpath,OUTN_VMOM_GRADKE,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hv,tmin,tmax,avg_start_time);
  VMom_dhdt = rho0*do_avg(dirpath,OUTN_VMOM_DHDT,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hv,tmin,tmax,avg_start_time);
  VMom_hypervisc = rho0*do_avg(dirpath,OUTN_VMOM_A4,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hv,tmin,tmax,avg_start_time);
  VMom_quadBotDrag = rho0*do_avg(dirpath,OUTN_VMOM_CDBOT,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hv,tmin,tmax,avg_start_time);
  VMom_windStress = rho0*do_avg(dirpath,OUTN_VMOM_WIND,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hv,tmin,tmax,avg_start_time);
  VMom_relaxation = rho0*do_avg(dirpath,OUTN_VMOM_RELAX,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hv,tmin,tmax,avg_start_time);
  VMom_baroForcing = rho0*do_avg(dirpath,OUTN_VMOM_FBARO,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hv,tmin,tmax,avg_start_time);
  VMom_randomForcing = rho0*do_avg(dirpath,OUTN_VMOM_RAND,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hv,tmin,tmax,avg_start_time);
  VMom_diapycnal = rho0*do_avg(dirpath,OUTN_VMOM_WDIA,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg_hv,tmin,tmax,avg_start_time);
  
  
  %%% Composite diagnostics
  UMom_advection = UMom_PVadvection + UMom_KEgradient + UMom_dhdt; 
  UMom_total = UMom_windStress+UMom_Montgomery+UMom_advection+UMom_quadBotDrag+UMom_hypervisc+UMom_relaxation+UMom_randomForcing+UMom_baroForcing+UMom_diapycnal;
  VMom_advection = VMom_PVadvection + VMom_KEgradient + VMom_dhdt;
  VMom_total = VMom_windStress+VMom_Montgomery+VMom_advection+VMom_quadBotDrag+VMom_hypervisc+VMom_relaxation+VMom_randomForcing+VMom_baroForcing+VMom_diapycnal;

  %%% Calculate streamfunctions
  hu = do_avg(dirpath,OUTN_HU_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
  Psi = zeros(Nx+1,Ny+1);
  Psi(1:Nx,2:Ny+1) = -cumsum(sum(hu,3)*dy,2);
  Psi(Nx+1,:) = Psi(1,:);
  Psi_uc = zeros(Nx+1,Ny+1);
  Psi_uc(1:Nx,2:Ny+1) = -cumsum(sum(hu(:,:,uc_layidx:end),3)*dy,2);
  Psi_uc(Nx+1,:) = Psi_uc(1,:);

  %%% Select variable whose isopleths we will use to average the circulation
  switch (coord_name)
    case 'bathy'
      coord_q = 0*hhb;
      coord_q(:,1) = 0.5*(hhb(:,1)+hhb([Nx 1:Nx-1],1));
      coord_q(:,2:Ny) = 0.25*(hhb(:,1:Ny-1)+hhb([Nx 1:Nx-1],1:Ny-1)+hhb(:,2:Ny)+hhb([Nx 1:Nx-1],2:Ny));
    case 'psi'
      coord_q = Psi_uc(1:Nx,1:Ny)/1e6;
  end

  %%% Calculate vorticity budget quantities
  curl_windStress = calc_curl(UMom_windStress,VMom_windStress,dx,dy);
  curl_quadBotDrag = calc_curl(UMom_quadBotDrag,VMom_quadBotDrag,dx,dy);
  curl_Montgomery = calc_curl(UMom_Montgomery,VMom_Montgomery,dx,dy);
  curl_hypervisc = calc_curl(UMom_hypervisc,VMom_hypervisc,dx,dy);
  curl_advection = calc_curl(UMom_advection,VMom_advection,dx,dy);
  curl_baroForcing = calc_curl(UMom_baroForcing,VMom_baroForcing,dx,dy);
  curl_randomForcing = calc_curl(UMom_randomForcing,VMom_randomForcing,dx,dy);
  curl_total = calc_curl(UMom_total,VMom_total,dx,dy);

  %%% Interpolate to q-points
  etab_q = 0*hhb;
  etab_q(:,1) = 0.5*(hhb(:,1)+hhb([Nx 1:Nx-1],1));
  etab_q(:,2:Ny) = 0.25*(hhb(:,1:Ny-1)+hhb([Nx 1:Nx-1],1:Ny-1)+hhb(:,2:Ny)+hhb([Nx 1:Nx-1],2:Ny));

  %%% Compute circulation tendencies along isobaths
  circ_windStress = calc_iso_circ(curl_windStress,dx,dy,dd,coord_q);
  circ_quadBotDrag = calc_iso_circ(curl_quadBotDrag,dx,dy,dd,coord_q);
  circ_Montgomery = calc_iso_circ(curl_Montgomery,dx,dy,dd,coord_q);
  circ_hypervisc = calc_iso_circ(curl_hypervisc,dx,dy,dd,coord_q);
  circ_advection = calc_iso_circ(curl_advection,dx,dy,dd,coord_q);
  circ_randomForcing = calc_iso_circ(curl_randomForcing,dx,dy,dd,coord_q);
  circ_baroForcing = calc_iso_circ(curl_baroForcing,dx,dy,dd,coord_q);
  circ_total = calc_iso_circ(curl_total,dx,dy,dd,coord_q);

  %%% Calculate contour lengths
  cntrlen = calc_iso_int(ones(Nx,Ny),dx,dy,dd,coord_q);

  %%% Calculate mean latitude of each isobath
  y_avg = calc_iso_int (YY_q(1:Nx,1:Ny),dx,dy,dd,coord_q) ./ cntrlen;
  
  save(fullfile('products',[run_name,'_momBalance.mat']), ...
    'dd','coord_name','Nd','coord_q', 'cntrlen','Psi_uc','y_avg', ...
    'UMom_windStress','UMom_Montgomery','UMom_advection','UMom_quadBotDrag','UMom_hypervisc','UMom_relaxation','UMom_randomForcing','UMom_baroForcing','UMom_diapycnal','UMom_total', ...
    'VMom_windStress','VMom_Montgomery','VMom_advection','VMom_quadBotDrag','VMom_hypervisc','VMom_relaxation','VMom_randomForcing','VMom_baroForcing','VMom_diapycnal','VMom_total', ... 
    'curl_windStress','curl_Montgomery','curl_advection','curl_quadBotDrag','curl_hypervisc','curl_randomForcing','curl_baroForcing','curl_total', ...
    'circ_windStress','circ_Montgomery','circ_advection','circ_quadBotDrag','circ_hypervisc','circ_randomForcing','circ_baroForcing','circ_total');

end

