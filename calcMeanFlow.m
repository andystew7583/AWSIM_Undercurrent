%%%
%%% calcMeanFlow.m
%%%
%%% Calculates diagnostics from time-averaged mean flow for one of our
%%% undercurrent runs.
%%%
function calcMeanFlow (local_home_dir,run_name,tmax,uc_layidx)

  %%% Define coordinate system for analysis
  coord_name = 'bathy';
  % coord_name = 'psi';

  %%% Load experiment parameters
  dirpath = fullfile(local_home_dir,run_name);
  loadParams;

  %%% Define coordinate grid
  switch (coord_name)
    case 'bathy'
      dd = [110:10:3500]';   
    case 'psi'
      dd = [0.01:0.01:0.5 0.55:0.05:5];
  end
  Nd = length(dd);

  %%% Average diagnostics
  tmax = tmax + 0.05*t1year;
  tmin = tmax - 5*t1year;
  u = do_avg(dirpath,OUTN_U_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
  v = do_avg(dirpath,OUTN_V_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
  h = do_avg(dirpath,OUTN_H_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
  hu = do_avg(dirpath,OUTN_HU_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
  hv = do_avg(dirpath,OUTN_HV_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);

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

  %%% Calculate curls of velocity and transport
  curl_u = calc_curl(u,v,dx,dy);
  curl_hu = calc_curl(hu,v,dx,dy);

  %%% Calculate contour lengths
  cntrlen = calc_iso_int (ones(Nx,Ny),dx,dy,dd,coord_q);

  %%% Calculate along-isobath circulation
  circ_u = calc_iso_circ(curl_u,dx,dy,dd,coord_q) ./ repmat(cntrlen,[1 Nlay]);
  circ_hu = calc_iso_circ(curl_hu,dx,dy,dd,coord_q) ./ repmat(cntrlen,[1 Nlay]);

  %%% Calculate mean latitude of each isobath
  y_avg = calc_iso_int (YY_q(1:Nx,1:Ny),dx,dy,dd,coord_q) ./ cntrlen;

  %%% Calculate along-isobath mean thickness
  h_avg = calc_iso_int (h_q,dx,dy,dd,coord_q) ./ repmat(cntrlen,[1 Nlay]);
  etab_avg = calc_iso_int (etab_q,dx,dy,dd,coord_q) ./ cntrlen;

  %%% Thickness-weighted average velocity
  circ_twa = circ_hu ./ h_avg;

  %%% Construct along-isobath mean layer elevations
  e_avg = zeros(Nd,Nlay+1);
  e_avg(:,Nlay+1) = etab_avg;
  for k=Nlay:-1:1
    e_avg(:,k) = e_avg(:,k+1) + h_avg(:,k);
  end

  save(fullfile('products',[run_name,'_meanFlow.mat']), ...
    'u','v','h','hu','hv','e','Psi_uc', ...
    'dd','Nd','coord_q', 'cntrlen', ...
    'curl_u','curl_hu', ...
    'circ_u','circ_hu', ...
    'y_avg','h_avg','etab_avg','e_avg', ...
    'circ_twa');

end

