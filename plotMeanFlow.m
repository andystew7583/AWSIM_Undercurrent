%%%
%%% plotMeanFlow.m
%%%
%%% Plots diagnostics from time-averaged mean flow for our
%%% undercurrent runs.
%%%

%%% Select experiment
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';
% run_name = 'undercurrent_RF_F0.05_Nc4_Ht200';
% run_name = 'undercurrent_WF_tau0.05_Nc4_localwind';
% run_name = 'undercurrent_WF_tau0.05_Nc4_Ht200';
% run_name = 'test_undercurrent_LWF_tau0.05_Nc4_Hbbl50_Nlay7';
% run_name = 'test_undercurrent_LWF_tau0.05_Nc4_Hbbl50_Nlay4';
% run_name = 'test_undercurrent_LWF_tau0.05_Nc4_Hbbl50_Nlay4_Htop250';
% run_name = 'test_undercurrent_LWF_tau0.05_Nc4_Hbbl50_Nlay4_Htop250_Huc400';
% run_name = 'test_undercurrent_LWF_tau0.05_Nc4_Hbbl50_Nlay7_Fbaro0.05';

% run_name = 'test_undercurrent_RF_F0.05_Nc4_Hbbl50';
% uc_layidx = 2; %%% Undercurrent layer
% tmin = 5.05*t1year;
% tmax = 15.05*t1year;

% run_name = 'test_undercurrent_RF_F0.05_Nc4_weakshelfforcing';
% uc_layidx = 2; %%% Undercurrent layer
% tmin = 5.05*t1year;
% tmax = 9.05*t1year;

% run_name = 'test_undercurrent_RF_F0.5_Nc4';
% uc_layidx = 2; %%% Undercurrent layer
% tmin = 5.05*t1year;
% tmax = 10.05*t1year;

run_name = 'test_undercurrent_RF_F0.5_Nc4_baroclinic';
uc_layidx = 2; %%% Undercurrent layer
tmin = 5.05*t1year;
tmax = 10.05*t1year;

% run_name = 'test_undercurrent_LWF_tau0.05_Nc4_Hbbl50_Nlay7';
% uc_layidx = 3; %%% Undercurrent layer
% tmin = 5.05*t1year;
% tmax = 10.05*t1year;

% run_name = 'test_undercurrent_LWF_tau0.05_Nc4_Hbbl50_Nlay7_Fbaro0.05';
% uc_layidx = 3; %%% Undercurrent layer
% tmin = 5.05*t1year;
% tmax = 10.05*t1year;

%%% Define coordinate system for analysis
% coord_name = 'bathy';
coord_name = 'psi';

%%% Load experiment parameters
dirpath = fullfile(local_home_dir,run_name);
loadParams;

%%% For writing figures
write_figs = true;
figdir = fullfile('plots',run_name);
mkdir(figdir);

%%% Define coordinate grid
switch (coord_name)
  case 'bathy'
    dd = [110:10:3500]';   
  case 'psi'
    dd = [0.01:0.01:0.5 0.55:0.05:5];
end
Nd = length(dd);

%%% Average diagnostics  
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

%%% Create grid for along-isobath averaged contour plots
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

%%% Create grid for slice contour plots
eps = 1; %%% Small distance by which to perturb plotting points away from actual isopycnal elevations
slice_idx = Nx/2;
% slice_idx = 3*Nx/8;
ZZ_slice = zeros(Ny,2*Nlay);
UU_slice = zeros(Ny,2*Nlay);
for k=1:Nlay
  ZZ_slice(:,2*k-1) = squeeze(e(slice_idx,:,k))-eps;
  ZZ_slice(:,2*k) = squeeze(e(slice_idx,:,k+1))+eps;
  UU_slice(:,2*k-1) = squeeze(u(slice_idx,:,k));
  UU_slice(:,2*k) = squeeze(u(slice_idx,:,k));
end
YY_slice = repmat(yy_h',[1 2*Nlay]);

%%% Plotting options
linewidth = 1.5;
axlim_zon = 0.2;
axlim_iso = 0.15;
ylim_iso = [y_avg(1)/1000 100];
ylim_zon = [0 350];
switch (coord_name)
  case 'bathy'
    d_levs = [110 150 200 400 1000 1500 2000 2500 3000 3500];
  case 'psi'
    d_levs = [0.1 0.5 1 3 5];
end
d_idx = zeros(1,length(d_levs));
d_labels = cell(1,length(d_levs));
for m=1:length(d_levs)
  d_idx(m) = find(abs(dd-d_levs(m))<1e-15);
  d_labels{m} = num2str(d_levs(m));
end
fontsize = 14;
axpos = [ 0.1300    0.1100    0.7750    0.73];
framepos = [527   853   560   420];
legstr = cell(1,Nlay);
for k=1:Nlay
  legstr{k} = ['layer ',num2str(k)];
end


figure(2);
clf;
set(gcf,'Position',framepos);
axes('Position',axpos);
plot(y_avg/1000,circ_u);
legend(legstr,'Location','NorthWest');
xlabel('Mean isobath offshore distance (km)');
set(gca,'XDir','reverse');
ylabel('m/s');
set(gca,'YLim',[-axlim_iso axlim_iso]);
set(gca,'XLim',ylim_iso);
set(gca,'FontSize',fontsize);
box off;
ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
plot(ax2,y_avg/1000,0*dd,'k--');
set(ax2,'YAxisLocation','Right');
set(ax2,'XAxisLocation','Top');
set(ax2,'XDir','reverse');
set(ax2,'Color','None');
set(ax2,'YTick',[]);
set(ax2,'XTick',y_avg(d_idx)/1000);
set(ax2,'XTickLabel',d_labels);
set(ax2,'YLim',[-axlim_iso axlim_iso]);
set(ax2,'XLim',ylim_iso);
set(ax2,'FontSize',fontsize);
xlabel(ax2,'Isobath (m)');
box off;
title(['Thickness-weighted average along-isobath mean flow']);
if (write_figs)
  print('-dpng','-r150',fullfile(figdir,'AlongIsobathMeanFlow.png'));
end

%%% Redefine axis position
axpos = [ 0.1300    0.1100    0.7750    0.8];


figure(4);
clf;
set(gcf,'Position',[359         175        1228         976]);
axes('Position',axpos);
pcolor(YY_avg/1000,-ZZ_avg,UU_avg);
shading interp
hold on;
plot(y_avg/1000,-e_avg,'Color',[.7 .7 .7]);
hold off;
colormap redblue(40)
colorbar;
caxis([-.1 .1])
xlabel('Mean isobath offshore distance (km)');
ylabel('Depth (m)');
set(gca,'XLim',[y_avg(1)/1000 100]);
set(gca,'YLim',[0 3000]);
set(gca,'XDir','reverse');
set(gca,'YDir','reverse');
set(gca,'FontSize',fontsize);
title('Thickness-weighted average along-isobath mean flow');
if (write_figs)
  print('-dpng','-r150',fullfile(figdir,'AlongIsobathFlowStructure.png'));
end

figure(5)
clf;
set(gcf,'Position',[874         713        1142         423]);
axes('Position',axpos);
pcolor(XX_u/1000,YY_u/1000,u(:,:,uc_layidx))
shading interp
hold on;
[C,h] = contour(XX_h/1000,YY_h/1000,-hhb,[200 1000 2000 3000 4000],'EdgeColor','k');
clabel(C,h);
hold off
colorbar
colormap redblue(40);
caxis([-axlim_zon axlim_zon]);
set(gca,'FontSize',fontsize);
xlabel('Alongshore distance (km)');
ylabel('Offshore distance (km)');
title(['Time-mean flow in layer ',num2str(uc_layidx)]);
if (write_figs)
  print('-dpng','-r150',fullfile(figdir,['AlongCoastFlow_layer',num2str(uc_layidx),'.png']));
end

figure(6);
clf;
set(gcf,'Position',[359         175        1228         976]);
axes('Position',axpos);
pcolor(YY_slice/1000,-ZZ_slice,UU_slice);
shading interp
hold on;
plot(yy_h/1000,-squeeze(e(slice_idx,:,:)),'Color',[.7 .7 .7]);
hold off;
colormap redblue(40)
colorbar;
caxis([-axlim_zon axlim_zon]);
xlabel('Offshore distance (km)');
ylabel('Depth (m)');
set(gca,'XDir','reverse');
set(gca,'YDir','reverse');
set(gca,'XLim',ylim_zon);
set(gca,'FontSize',fontsize);
title(['Time-mean along-coast flow at x=',num2str(round(xx_h(slice_idx)/1000)),' km']);
if (write_figs)
  print('-dpng','-r150',fullfile(figdir,['AlongCoastFlow_x=',num2str(round(xx_h(slice_idx)/1000)),'.png']));
end

figure(7);
pcolor(XX_q/1000,YY_q/1000,Psi/1e6);
shading interp;
colormap redblue;
colorbar;

figure(8);
pcolor(XX_q/1000,YY_q/1000,Psi_uc/1e6);
shading interp;
colormap redblue;
colorbar;