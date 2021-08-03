%%%
%%% plotMeanFlow.m
%%%
%%% Plots diagnostics from time-averaged mean flow for our
%%% undercurrent runs.
%%%

%%% Select experiment
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';
% run_name = 'undercurrent_RF_F0.05_Nc4_Ht200';
run_name = 'undercurrent_WF_tau0.05_Nc4_localwind';
% run_name = 'undercurrent_WF_tau0.05_Nc4_Ht200';

%%% Load experiment parameters
dirpath = fullfile(local_home_dir,run_name);
loadParams;
tmin = 5.05*t1year;
tmax = 10.05*t1year;

%%% Define isobath grid
dd = [100:10:4000]';
Nd = length(dd);

%%% Average diagnostics  
u = do_avg(dirpath,OUTN_U_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
v = do_avg(dirpath,OUTN_V_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
h = do_avg(dirpath,OUTN_H_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
hu = do_avg(dirpath,OUTN_HU_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);
hv = do_avg(dirpath,OUTN_HV_AVG,Nx,Ny,Nlay,n0_avg,N_avg,dt_avg,tmin,tmax,startTime);

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

curl_u = calc_curl(u,v,dx,dy);
curl_hu = calc_curl(hu,v,dx,dy);

%%% Calculate contour lengths
cntrlen = calc_iso_int (ones(Nx,Ny),dx,dy,dd,etab_q);

%%% Calculate along-isobath circulation
circ_u = calc_iso_circ(curl_u,dx,dy,dd,etab_q) ./ repmat(cntrlen,[1 Nlay]);
circ_hu = calc_iso_circ(curl_hu,dx,dy,dd,etab_q) ./ repmat(cntrlen,[1 Nlay]);

%%% Calculate mean latitude of each isobath
y_avg = calc_iso_int (YY_q(1:Nx,1:Ny),dx,dy,dd,etab_q) ./ cntrlen;

%%% Calculate along-isobath mean thickness
h_avg = calc_iso_int (h_q,dx,dy,dd,etab_q) ./ repmat(cntrlen,[1 Nlay]);

%%% Thickness-weighted average velocity
circ_twa = circ_hu ./ h_avg;

%%% Construct along-isobath mean layer elevations
e_avg = zeros(Nd,Nlay+1);
e_avg(:,Nlay+1) = -dd;
for k=Nlay:-1:1
  e_avg(:,k) = e_avg(:,k+1) + h_avg(:,k);
end

%%% Create grid for contour plots
eps = 1; %%% Small distance by which to perturb plotting points away from actual isopycnal elevations
ZZ_avg = zeros(Nd,2*Nlay);
UU_avg = zeros(Nd,2*Nlay);
for k=1:Nlay
  ZZ_avg(:,2*k-1) = e_avg(:,k)-eps;
  ZZ_avg(:,2*k) = e_avg(:,k+1)+eps;
  UU_avg(:,2*k-1) = circ_twa(:,k);
  UU_avg(:,2*k) = circ_twa(:,k);
end
YY_avg = repmat(y_avg,[1 2*Nlay]);

figure(1);
plot(yy_h,squeeze(mean(u,1)));

figure(2);
plot(dd,circ_u);

figure(3)
plot(y_avg,e_avg);

figure(4);
pcolor(YY_avg,ZZ_avg,UU_avg);
shading interp
hold on;
plot(y_avg,e_avg,'Color',[.7 .7 .7]);
hold off;
colormap redblue
colorbar;
caxis([-.1 .1])
