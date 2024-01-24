%%%
%%% calcUCflowTimeSeries.m
%%%
%%% Calculates a time series of simulated undercurrent flow speed.
%%%

%%% Load parameters   
loadParams;
dirpath = fullfile(local_home_dir,run_name);

%%% Define coordinates for averaging
coord_q = 0*hhb;
coord_q(:,1) = 0.5*(hhb(:,1)+hhb([Nx 1:Nx-1],1));
coord_q(:,2:Ny) = 0.25*(hhb(:,1:Ny-1)+hhb([Nx 1:Nx-1],1:Ny-1)+hhb(:,2:Ny)+hhb([Nx 1:Nx-1],2:Ny));
dd = [110:10:3500]';   
uc_didx = find((dd>=120) & (dd<=600));
uc_layidx = 3;

%%% Calculate contour lengths
cntrlen = calc_iso_int (ones(Nx,Ny),dx,dy,dd,coord_q);

Mcnt = 1;

n0 = 0;
Nframes = endTime/dt_s;
startTime = 0;

%%% At each time iteration...
for n=n0:1:n0+Nframes-1   

  %%% Current simulation time    
  t = startTime + (n-n0)*dt_s;
  tt(Mcnt) = t;  

  disp(t/t1year)

  %%% Load model state
  layer = uc_layidx;
  data_file = fullfile(dirpath,[OUTN_U,num2str(layer-1),'_n=',num2str(n),'.dat']);
  uu = readOutputFile(data_file,Nx,Ny);     
  data_file = fullfile(dirpath,[OUTN_V,num2str(layer-1),'_n=',num2str(n),'.dat']);
  vv = readOutputFile(data_file,Nx,Ny);     
  data_file = fullfile(dirpath,[OUTN_H,num2str(layer-1),'_n=',num2str(n),'.dat']);
  h = readOutputFile(data_file,Nx,Ny);     
%     data_file = fullfile(dirpath,[OUTN_HU_AVG,num2str(layer-1),'_n=',num2str(n),'.dat']);
%     hu = readOutputFile(data_file,Nx,Ny);    
%     data_file = fullfile(dirpath,[OUTN_HV_AVG,num2str(layer-1),'_n=',num2str(n),'.dat']);
%     hv = readOutputFile(data_file,Nx,Ny);   
%     data_file = fullfile(dirpath,[OUTN_H_AVG,num2str(layer-1),'_n=',num2str(n),'.dat']);
%     h = readOutputFile(data_file,Nx,Ny);     
  hu = 0.5*(h(1:Nx,1:Ny)+h([Nx 1:Nx-1],1:Ny)).*uu;
  hv = 0.5*(h(1:Nx,1:Ny)+h(1:Nx,[Ny 1:Ny-1])).*vv;

  if (isempty(h))      
      u_UC(Mcnt) = NaN;
      continue;
  end


  h_q = 0*h;
  h_q(:,1,:) = 0.5*(h(:,1,:)+h([Nx 1:Nx-1],1,:));
  h_q(:,2:Ny,:) = 0.25*(h(:,1:Ny-1,:)+h([Nx 1:Nx-1],1:Ny-1,:)+h(:,2:Ny,:)+h([Nx 1:Nx-1],2:Ny,:));


  %%% Calculate curls of velocity and transport
%     curl_u = calc_curl(uu,vv,dx,dy);
  curl_hu = calc_curl(hu,hv,dx,dy);    

  %%% Calculate along-isobath circulation
%     circ_u = calc_iso_circ(curl_u,dx,dy,dd,coord_q) ./ repmat(cntrlen,[1 Nlay]);
  circ_hu = calc_iso_circ(curl_hu,dx,dy,dd,coord_q) ./ repmat(cntrlen,[1 Nlay]);

  %%% Thickness-weighted averaged along-isobath flow
  h_avg = calc_iso_int (h_q,dx,dy,dd,coord_q) ./ repmat(cntrlen,[1 Nlay]);    
  u_twa = circ_hu ./ h_avg;

  %%% Average over isobaths to get undercurrent speed
  u_UC(Mcnt) = mean(u_twa(uc_didx));

  Mcnt = Mcnt + 1;

end

save([run_name,'_UCflow.mat'],'tt','u_UC');
