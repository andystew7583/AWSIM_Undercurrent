%%%
%%% calcMeanFlow.m
%%%
%%% Calculates diagnostics from time-averaged mean flow for one of our
%%% undercurrent runs.
%%%
function calcMeanFlow (local_home_dir,run_name,tmax,uc_layidx)

  %%% Define coordinate system for analysis
  coord_name = 'bathy';
%   coord_name = 'psi';

  %%% Load experiment parameters
  dirpath = fullfile(local_home_dir,run_name);
  loadParams;

  %%% Define coordinate grid
  switch (coord_name)
    case 'bathy'
      dd = [110:10:3500]';   
    case 'psi'
      dd = [0.001:0.001:0.01 0.02:0.01:0.5 0.55:0.05:5];
  end
  Nd = length(dd);

  %%% Average diagnostics
%   tmax = tmax + 0.05*t1year;
  tmax = tmax - 2*t1year+ 0.05*t1year;
%   tmin = tmax - 10*t1year;
  tmin = tmax - 8*t1year;

%   avg_iter_start = n0_avg;
%   avg_num_iters = N_avg;
%   avg_start_time = startTime;
  avg_iter_start = 0;                    
  avg_num_iters = N_avg+n0_avg;
  avg_start_time = 0;
  u = do_avg(dirpath,OUTN_U_AVG,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg,tmin,tmax,avg_start_time);
  v = do_avg(dirpath,OUTN_V_AVG,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg,tmin,tmax,avg_start_time);
  h = do_avg(dirpath,OUTN_H_AVG,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg,tmin,tmax,avg_start_time);
  M = do_avg(dirpath,OUTN_M_AVG,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg,tmin,tmax,avg_start_time);
  hu = do_avg(dirpath,OUTN_HU_AVG,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg,tmin,tmax,avg_start_time);
  hv = do_avg(dirpath,OUTN_HV_AVG,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg,tmin,tmax,avg_start_time);
  huu = do_avg(dirpath,OUTN_HUU_AVG,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg,tmin,tmax,avg_start_time);
  hvv = do_avg(dirpath,OUTN_HVV_AVG,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg,tmin,tmax,avg_start_time);  
  huv = do_avg(dirpath,OUTN_HUV_AVG,Nx,Ny,Nlay,avg_iter_start,avg_num_iters,dt_avg,tmin,tmax,avg_start_time);  
  

  %%% Calculate streamfunctions
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

  %%% Mean layer interfaces
  e = zeros(Nx,Ny,Nlay+1);
  e(:,:,end) = hhb;
  for k=Nlay:-1:1
    e(:,:,k) = e(:,:,k+1) + h(:,:,k);
  end

  %%% Interpolate to q-points
  etab_q = 0*hhb;
  etab_q(:,1) = 0.5*(hhb(:,1)+hhb([Nx 1:Nx-1],1));
  etab_q(:,2:Ny) = 0.25*(hhb(:,1:Ny-1)+hhb([Nx 1:Nx-1],1:Ny-1)+hhb(:,2:Ny)+hhb([Nx 1:Nx-1],2:Ny));
  h_q = 0*h;
  h_q(:,1,:) = 0.5*(h(:,1,:)+h([Nx 1:Nx-1],1,:));
  h_q(:,2:Ny,:) = 0.25*(h(:,1:Ny-1,:)+h([Nx 1:Nx-1],1:Ny-1,:)+h(:,2:Ny,:)+h([Nx 1:Nx-1],2:Ny,:));
  hu_q = 0*hu;
  hu_q(:,2:Ny,:) = 0.5*(hu(:,1:Ny-1,:)+hu(:,2:Ny,:));
  hu_q(:,1,:) = hu(:,1,:);
  huu_q = 0*huu;
  huu_q(:,2:Ny,:) = 0.5*(huu(:,1:Ny-1,:)+huu(:,2:Ny,:));
  huu_q(:,1,:) = huu(:,1,:);
  hv_q = 0.5*(hv(1:Nx,:,:)+hv([Nx 1:Nx-1],:,:));
  hvv_q = 0.5*(hvv(1:Nx,:,:)+hvv([Nx 1:Nx-1],:,:));

  %%% Calculate curls of velocity and transport
  curl_u = calc_curl(u,v,dx,dy);
  curl_hu = calc_curl(hu,hv,dx,dy);
 
  %%% Calculate contour lengths
  cntrlen = calc_iso_int (ones(Nx,Ny),dx,dy,dd,coord_q);

  %%% Calculate along-isobath circulation
  circ_u = calc_iso_circ(curl_u,dx,dy,dd,coord_q) ./ repmat(cntrlen,[1 Nlay]);
  circ_hu = calc_iso_circ(curl_hu,dx,dy,dd,coord_q) ./ repmat(cntrlen,[1 Nlay]);
  


  %%% Calculate mean latitude of each isobath
  y_avg = calc_iso_int (YY_q(1:Nx,1:Ny),dx,dy,dd,coord_q) ./ cntrlen;

  %%% Calculate along-isobath mean quantities
  h_avg = calc_iso_int (h_q,dx,dy,dd,coord_q) ./ repmat(cntrlen,[1 Nlay]);
  etab_avg = calc_iso_int (etab_q,dx,dy,dd,coord_q) ./ cntrlen;
  hu_avg = calc_iso_int (hu_q,dx,dy,dd,coord_q) ./ repmat(cntrlen,[1 Nlay]);
  huu_avg = calc_iso_int (huu_q,dx,dy,dd,coord_q) ./ repmat(cntrlen,[1 Nlay]);
  hv_avg = calc_iso_int (hv_q,dx,dy,dd,coord_q) ./ repmat(cntrlen,[1 Nlay]);
  hvv_avg = calc_iso_int (hvv_q,dx,dy,dd,coord_q) ./ repmat(cntrlen,[1 Nlay]);
  
  %%% Construct along-isobath mean layer elevations
  e_avg = zeros(Nd,Nlay+1);
  e_avg(:,Nlay+1) = etab_avg;
  for k=Nlay:-1:1
    e_avg(:,k) = e_avg(:,k+1) + h_avg(:,k);
  end
  
  %%% Thickness-weighted average velocity
  circ_twa = circ_hu ./ h_avg;
  
  %%% Calculate EKE
  KE = huu_avg + hvv_avg;
  MKE = (hu_avg.^2 + hv_avg.^2)./h_avg;
  EKE = KE - MKE;

  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Calculate advective forcing %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  %%% Calculate mean and eddy momentum fluxes
  huu_mean = hu.^2./h;
  hvv_mean = hv.^2./h;
  huv_mean = hu.*hv./h;
  huu_eddy = huu - huu_mean;
  hvv_eddy = hvv - hvv_mean;
  huv_eddy = huv - huv_mean;
  
  %%% Calculate advective tendencies
  UMom_totalMomFlux = - (huu([2:Nx 1],:,:) - huu([Nx 1:Nx-1],:,:)) / (2*dx) ...
                      - (huv(:,[2:Ny 1],:) - huv(:,1:Ny,:)) / dy;
  VMom_totalMomFlux = - (huv([2:Nx 1],:,:) - huv(1:Nx,:,:)) / dx ...
                      - (hvv(:,[2:Ny 1],:) - hvv(:,[Ny 1:Ny-1],:)) / (2*dy); 
  UMom_meanMomFlux = - (huu_mean([2:Nx 1],:,:) - huu_mean([Nx 1:Nx-1],:,:)) / (2*dx) ...
                      - (huv_mean(:,[2:Ny 1],:) - huv_mean(:,1:Ny,:)) / dy;
  VMom_meanMomFlux = - (huv_mean([2:Nx 1],:,:) - huv_mean(1:Nx,:,:)) / dx ...
                      - (hvv_mean(:,[2:Ny 1],:) - hvv_mean(:,[Ny 1:Ny-1],:)) / (2*dy); 
  UMom_eddyMomFlux = - (huu_eddy([2:Nx 1],:,:) - huu_eddy([Nx 1:Nx-1],:,:)) / (2*dx) ...
                      - (huv_eddy(:,[2:Ny 1],:) - huv_eddy(:,1:Ny,:)) / dy;
  VMom_eddyMomFlux = - (huv_eddy([2:Nx 1],:,:) - huv_eddy(1:Nx,:,:)) / dx ...
                      - (hvv_eddy(:,[2:Ny 1],:) - hvv_eddy(:,[Ny 1:Ny-1],:)) / (2*dy); 
    
  %%% Calculate curl of advective forcing
  curl_totalMomFlux = calc_curl(UMom_totalMomFlux,VMom_totalMomFlux,dx,dy);
  curl_meanMomFlux = calc_curl(UMom_meanMomFlux,VMom_meanMomFlux,dx,dy);
  curl_eddyMomFlux = calc_curl(UMom_eddyMomFlux,VMom_eddyMomFlux,dx,dy);

  %%% Compute advective circulation tendencies along isobaths
  circ_totalMomFlux = calc_iso_circ(curl_totalMomFlux,dx,dy,dd,coord_q);
  circ_meanMomFlux = calc_iso_circ(curl_meanMomFlux,dx,dy,dd,coord_q);
  circ_eddyMomFlux = calc_iso_circ(curl_eddyMomFlux,dx,dy,dd,coord_q);
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%% Calculate pressure forcing %%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %%% TODO decompose h grad M
  
  

  save(fullfile('products',[run_name,'_meanFlow.mat']), ...
    'u','v','h','hu','hv','huu','hvv','huv','e','Psi','Psi_uc', ...
    'coord_name','dd','Nd','coord_q', 'cntrlen', ...
    'curl_u','curl_hu', ...
    'circ_u','circ_hu', ...
    'y_avg','h_avg','etab_avg','e_avg', ...
    'hu_avg','huu_avg','hv_avg','hvv_avg', ...
    'KE','MKE','EKE', ...
    'circ_twa', ...
    'curl_totalMomFlux','curl_meanMomFlux','curl_eddyMomFlux', ...
    'circ_totalMomFlux','circ_meanMomFlux','circ_eddyMomFlux');

end

