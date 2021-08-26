%%%
%%% plotMomBudget.m
%%%
%%% Plots diagnostics from time-averaged momentum budget for our
%%% undercurrent runs.
%%%

%%% Select experiment
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';
% run_name = 'undercurrent_RF_F0.05_Nc4_Ht200';
% run_name = 'undercurrent_WF_tau0.05_Nc4_localwind';
% run_name = 'undercurrent_WF_tau0.05_Nc4_Ht200';
% run_name = 'test_undercurrent_LWF_tau0.05_Nc4_Hbbl50_Nlay4_Htop250';
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

%%% Load experiment parameters
dirpath = fullfile(local_home_dir,run_name);
loadParams;
rho0 = 1000;

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

%%% Average u-momentum diagnostics  
UMom_PVadvection = rho0*do_avg(dirpath,OUTN_UMOM_Q,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_Montgomery = rho0*do_avg(dirpath,OUTN_UMOM_GRADM,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_KEgradient = rho0*do_avg(dirpath,OUTN_UMOM_GRADKE,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_dhdt = rho0*do_avg(dirpath,OUTN_UMOM_DHDT,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_hypervisc = rho0*do_avg(dirpath,OUTN_UMOM_A4,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_quadBotDrag = rho0*do_avg(dirpath,OUTN_UMOM_CDBOT,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_windStress = rho0*do_avg(dirpath,OUTN_UMOM_WIND,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_relaxation = rho0*do_avg(dirpath,OUTN_UMOM_RELAX,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_baroForcing = rho0*do_avg(dirpath,OUTN_UMOM_FBARO,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_randomForcing = rho0*do_avg(dirpath,OUTN_UMOM_RAND,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_diapycnal = rho0*do_avg(dirpath,OUTN_UMOM_WDIA,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
% UMom_buoyForce = rho0*do_avg(dirpath,OUTN_UMOM_BUOY,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
% UMom_linBotDrag = rho0*do_avg(dirpath,OUTN_UMOM_RDRAG,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
% UMom_linSurfDrag = rho0*do_avg(dirpath,OUTN_UMOM_RSURF,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
% UMom_quadSurfDrag = rho0*do_avg(dirpath,OUTN_UMOM_CDSURF,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_advection = UMom_PVadvection + UMom_KEgradient + UMom_dhdt; 
UMom_total = UMom_windStress+UMom_Montgomery+UMom_advection+UMom_quadBotDrag+UMom_hypervisc+UMom_relaxation+UMom_randomForcing+UMom_baroForcing+UMom_diapycnal;

%%% Average v-momentum diagnostics  
VMom_PVadvection = rho0*do_avg(dirpath,OUTN_VMOM_Q,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_Montgomery = rho0*do_avg(dirpath,OUTN_VMOM_GRADM,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_KEgradient = rho0*do_avg(dirpath,OUTN_VMOM_GRADKE,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_dhdt = rho0*do_avg(dirpath,OUTN_VMOM_DHDT,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_hypervisc = rho0*do_avg(dirpath,OUTN_VMOM_A4,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_quadBotDrag = rho0*do_avg(dirpath,OUTN_VMOM_CDBOT,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_windStress = rho0*do_avg(dirpath,OUTN_VMOM_WIND,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_relaxation = rho0*do_avg(dirpath,OUTN_VMOM_RELAX,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_baroForcing = rho0*do_avg(dirpath,OUTN_VMOM_FBARO,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_randomForcing = rho0*do_avg(dirpath,OUTN_VMOM_RAND,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_diapycnal = rho0*do_avg(dirpath,OUTN_VMOM_WDIA,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
% VMom_buoyForce = rho0*do_avg(dirpath,OUTN_VMOM_BUOY,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
% VMom_linBotDrag = rho0*do_avg(dirpath,OUTN_VMOM_RDRAG,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
% VMom_linSurfDrag = rho0*do_avg(dirpath,OUTN_VMOM_RSURF,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
% VMom_quadSurfDrag = rho0*do_avg(dirpath,OUTN_VMOM_CDSURF,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
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





%%% Plotting options
linewidth = 1.5;
axlim_zon = 7e-2;
axlim_iso = 5e-2;
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

figure(10);
clf;
set(gcf,'Position',framepos);
axes('Position',axpos);
plot(yy_h/1000,mean(sum(UMom_windStress,3),1),'LineWidth',linewidth);
hold on;
plot(yy_h/1000,mean(sum(UMom_Montgomery,3),1),'LineWidth',linewidth);
plot(yy_h/1000,mean(sum(UMom_advection,3),1),'LineWidth',linewidth);
plot(yy_h/1000,mean(sum(UMom_quadBotDrag,3),1),'LineWidth',linewidth);
plot(yy_h/1000,mean(sum(UMom_hypervisc,3),1),'LineWidth',linewidth);
plot(yy_h/1000,mean(sum(UMom_randomForcing,3),1),'LineWidth',linewidth);
plot(yy_h/1000,mean(sum(UMom_baroForcing,3),1),'LineWidth',linewidth);
plot(yy_h/1000,mean(sum(UMom_total,3),1),'LineWidth',linewidth,'Color',[.7 .7 .7]);
hold off
legend('Wind','PGF','Advection','Friction','Viscosity','Random forcing','Barotropic forcing','Total','Location','NorthWest');
xlabel('Offshore distance (km)');
set(gca,'XDir','reverse');
ylabel('N/m^2');
set(gca,'YLim',[-axlim_zon axlim_zon]);
set(gca,'XLim',ylim_zon);
set(gca,'FontSize',fontsize);
title(['Along-shore momentum budget over full ocean depth']);
if (write_figs)
  print('-dpng','-r150',fullfile(figdir,'LinearMomBalance_FullDepth.png'));
end

figure(11);
clf;
set(gcf,'Position',framepos);
axes('Position',axpos);
plot(yy_h/1000,mean(sum(UMom_windStress(:,:,1:uc_layidx-1),3),1),'LineWidth',linewidth);
hold on;
plot(yy_h/1000,mean(sum(UMom_Montgomery(:,:,1:uc_layidx-1),3),1),'LineWidth',linewidth);
plot(yy_h/1000,mean(sum(UMom_advection(:,:,1:uc_layidx-1),3),1),'LineWidth',linewidth);
plot(yy_h/1000,mean(sum(UMom_quadBotDrag(:,:,1:uc_layidx-1),3),1),'LineWidth',linewidth);
plot(yy_h/1000,mean(sum(UMom_hypervisc(:,:,1:uc_layidx-1),3),1),'LineWidth',linewidth);
plot(yy_h/1000,mean(sum(UMom_randomForcing(:,:,1:uc_layidx-1),3),1),'LineWidth',linewidth);
plot(yy_h/1000,mean(sum(UMom_baroForcing(:,:,1:uc_layidx-1),3),1),'LineWidth',linewidth);
plot(yy_h/1000,mean(sum(UMom_total(:,:,1:uc_layidx-1),3),1),'LineWidth',linewidth,'Color',[.7 .7 .7]);
hold off
legend('Wind','PGF','Advection','Friction','Viscosity','Random forcing','Barotropic forcing','Total','Location','NorthWest');
xlabel('Offshore distance (km)');
set(gca,'XDir','reverse');
ylabel('N/m^2');
set(gca,'YLim',[-axlim_zon axlim_zon]);
set(gca,'XLim',ylim_zon);
set(gca,'FontSize',fontsize);
title(['Along-shore momentum budget, layers 1-',num2str(uc_layidx-1)]);
if (write_figs)
  print('-dpng','-r150',fullfile(figdir,'LinearMomBalance_Surface.png'));
end

figure(12);
clf;
set(gcf,'Position',framepos);
axes('Position',axpos);
plot(yy_h/1000,mean(sum(UMom_windStress(:,:,uc_layidx:Nlay),3),1),'LineWidth',linewidth);
hold on;
plot(yy_h/1000,mean(sum(UMom_Montgomery(:,:,uc_layidx:Nlay),3),1),'LineWidth',linewidth);
plot(yy_h/1000,mean(sum(UMom_advection(:,:,uc_layidx:Nlay),3),1),'LineWidth',linewidth);
plot(yy_h/1000,mean(sum(UMom_quadBotDrag(:,:,uc_layidx:Nlay),3),1),'LineWidth',linewidth);
plot(yy_h/1000,mean(sum(UMom_hypervisc(:,:,uc_layidx:Nlay),3),1),'LineWidth',linewidth);
plot(yy_h/1000,mean(sum(UMom_randomForcing(:,:,uc_layidx:Nlay),3),1),'LineWidth',linewidth);
plot(yy_h/1000,mean(sum(UMom_baroForcing(:,:,uc_layidx:Nlay),3),1),'LineWidth',linewidth);
plot(yy_h/1000,mean(sum(UMom_total(:,:,uc_layidx:Nlay),3),1),'LineWidth',linewidth,'Color',[.7 .7 .7]);
hold off
legend('Wind','PGF','Advection','Friction','Viscosity','Random forcing','Barotropic forcing','Total','Location','NorthWest');
xlabel('Offshore distance (km)');
set(gca,'XDir','reverse');
ylabel('N/m^2');
set(gca,'YLim',[-axlim_zon axlim_zon]);
set(gca,'XLim',ylim_zon);
set(gca,'FontSize',fontsize);
title(['Along-shore momentum budget, layers ',num2str(uc_layidx),'-',num2str(Nlay)]);
if (write_figs)
  print('-dpng','-r150',fullfile(figdir,'LinearMomBalance_Undercurrent.png'));
end

figure(13);
clf;
set(gcf,'Position',framepos);
axes('Position',axpos);
plot(y_avg/1000,sum(circ_windStress,2)./cntrlen,'LineWidth',linewidth);
hold on;
plot(y_avg/1000,sum(circ_Montgomery,2)./cntrlen,'LineWidth',linewidth);
plot(y_avg/1000,sum(circ_advection,2)./cntrlen,'LineWidth',linewidth);
plot(y_avg/1000,sum(circ_quadBotDrag,2)./cntrlen,'LineWidth',linewidth);
plot(y_avg/1000,sum(circ_hypervisc,2)./cntrlen,'LineWidth',linewidth);
plot(y_avg/1000,sum(circ_randomForcing,2)./cntrlen,'LineWidth',linewidth);
plot(y_avg/1000,sum(circ_baroForcing,2)./cntrlen,'LineWidth',linewidth);
plot(y_avg/1000,sum(circ_total,2)./cntrlen,'LineWidth',linewidth,'Color',[.7 .7 .7]);
hold off;
legend('Wind','PGF','Advection','Friction','Viscosity','Random forcing','Barotropic forcing','Total','Location','NorthWest');
xlabel('Mean isobath offshore distance (km)');
set(gca,'XDir','reverse');
ylabel('N/m^2');
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
title(['Along-isobath momentum budget over full ocean depth']);
if (write_figs)
  print('-dpng','-r150',fullfile(figdir,'IsobathMomBalance_FullDepth.png'));
end

figure(14);
clf;
set(gcf,'Position',framepos);
axes('Position',axpos);
plot(y_avg/1000,sum(circ_windStress(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth);
hold on;
plot(y_avg/1000,sum(circ_Montgomery(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth);
plot(y_avg/1000,sum(circ_advection(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth);
plot(y_avg/1000,sum(circ_quadBotDrag(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth);
plot(y_avg/1000,sum(circ_hypervisc(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth);
plot(y_avg/1000,sum(circ_randomForcing(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth);
plot(y_avg/1000,sum(circ_baroForcing(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth);
plot(y_avg/1000,sum(circ_total(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth,'Color',[.7 .7 .7]);
hold off;
legend('Wind','PGF','Advection','Friction','Viscosity','Random forcing','Barotropic forcing','Total','Location','NorthWest');
xlabel('Mean isobath offshore distance (km)');
set(gca,'XDir','reverse');
ylabel('N/m^2');
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
title(['Along-isobath momentum budget, layers 1-',num2str(uc_layidx-1)]);
if (write_figs)
  print('-dpng','-r150',fullfile(figdir,'IsobathMomBalance_Surface.png'));
end

figure(15);
clf;
set(gcf,'Position',framepos);
axes('Position',axpos);
plot(y_avg/1000,sum(circ_windStress(:,uc_layidx:Nlay),2)./cntrlen,'LineWidth',linewidth);
hold on;
plot(y_avg/1000,sum(circ_Montgomery(:,uc_layidx:Nlay),2)./cntrlen,'LineWidth',linewidth);
plot(y_avg/1000,sum(circ_advection(:,uc_layidx:Nlay),2)./cntrlen,'LineWidth',linewidth);
plot(y_avg/1000,sum(circ_quadBotDrag(:,uc_layidx:Nlay),2)./cntrlen,'LineWidth',linewidth);
plot(y_avg/1000,sum(circ_hypervisc(:,uc_layidx:Nlay),2)./cntrlen,'LineWidth',linewidth);
plot(y_avg/1000,sum(circ_randomForcing(:,uc_layidx:Nlay),2)./cntrlen,'LineWidth',linewidth);
plot(y_avg/1000,sum(circ_baroForcing(:,uc_layidx:Nlay),2)./cntrlen,'LineWidth',linewidth);
plot(y_avg/1000,sum(circ_total(:,uc_layidx:Nlay),2)./cntrlen,'LineWidth',linewidth,'Color',[.7 .7 .7]);
hold off;
legend('Wind','PGF','Advection','Friction','Viscosity','Random forcing','Barotropic forcing','Total','Location','NorthWest');
xlabel('Mean isobath offshore distance (km)');
set(gca,'XDir','reverse');
ylabel('N/m^2');
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
title(['Along-isobath momentum budget, layers ',num2str(uc_layidx),'-',num2str(Nlay)]);
if (write_figs)
  print('-dpng','-r150',fullfile(figdir,'IsobathMomBalance_Undercurrent.png'));
end