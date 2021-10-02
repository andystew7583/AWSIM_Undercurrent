%%%
%%% plotMeanFlow.m
%%%
%%% Plots diagnostics from time-averaged mean flow for our
%%% undercurrent runs.
%%%

%%% Location of runs on file system
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';

%%% Select model configuration and parameters
config = 'rand';
grid_size = 128; %%% Default 128
wind_stress = 0.05; %%% Default 0.05
rand_force = 0.75; %%% Default 0.75
num_canyons = 4; %%% Default 4
amp_canyons = 25; %%% Default 25
max_slope = 0.15; %%% Default 0.15
sb_width = 5; %%% Default 5
baro_force = 0; %%% Default 0
drag_coeff = 2; %%% Default 2

%%% Undercurrent layer
switch (config)
  case 'wind'
    uc_layidx = 3;
  case 'rand'
    uc_layidx = 2;
end

%%% Generate simulation name
run_name = constructRunName (config,grid_size,wind_stress, ...
          rand_force,num_canyons,amp_canyons,max_slope,sb_width,baro_force,drag_coeff);

%%% Load along-isobath mean flow diagnostics
loadParams;
rho0 = 1000;
load(fullfile('products',[run_name,'_meanFlow.mat']));

%%% For writing figures
write_figs = true;
figdir = fullfile('plots',run_name);
mkdir(figdir);

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

figure(9);
plot(y_avg/1000,EKE./h_avg);
xlabel('Mean isobath offshore distance (km)');
ylabel('EKE density (m^2/s^2)');
set(gca,'XDir','reverse');
set(gca,'XLim',ylim_iso);
set(gca,'FontSize',fontsize);

figure(10);
plot(y_avg/1000,KE./h_avg);
xlabel('Mean isobath offshore distance (km)');
ylabel('TKE density (m^2/s^2)');
set(gca,'XDir','reverse');
set(gca,'XLim',ylim_iso);
set(gca,'FontSize',fontsize);

figure(11);
plot(y_avg/1000,MKE./h_avg);
xlabel('Mean isobath offshore distance (km)');
ylabel('MKE density (m^2/s^2)');
set(gca,'XDir','reverse');
set(gca,'XLim',ylim_iso);
set(gca,'FontSize',fontsize);








%%% For momentum balance plots
axlim_iso = 5e-2;

figure(33);
clf;
set(gcf,'Position',framepos);
axes('Position',axpos);
plot(y_avg/1000,rho0*sum(circ_totalMomFlux,2)./cntrlen,'LineWidth',linewidth);
hold on;
plot(y_avg/1000,rho0*sum(circ_meanMomFlux,2)./cntrlen,'LineWidth',linewidth);
plot(y_avg/1000,rho0*sum(circ_eddyMomFlux,2)./cntrlen,'LineWidth',linewidth);
hold off;
legend('Total','Mean','Eddy','Location','NorthWest');
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
title(['Along-isobath advective momentum forcing over full ocean depth']);
if (write_figs)
  print('-dpng','-r150',fullfile(figdir,'IsobathAdvForcing_FullDepth.png'));
end

figure(34);
clf;
set(gcf,'Position',framepos);
axes('Position',axpos);
plot(y_avg/1000,rho0*sum(circ_totalMomFlux(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth);
hold on;
plot(y_avg/1000,rho0*sum(circ_meanMomFlux(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth);
plot(y_avg/1000,rho0*sum(circ_eddyMomFlux(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth);
hold off;
legend('Total','Mean','Eddy','Location','NorthWest');
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
title(['Along-isobath advective momentum forcing, layers 1-',num2str(uc_layidx-1)]);
if (write_figs)
  print('-dpng','-r150',fullfile(figdir,'IsobathAdvForcing_Surface.png'));
end

figure(35);
clf;
set(gcf,'Position',framepos);
axes('Position',axpos);
plot(y_avg/1000,rho0*sum(circ_totalMomFlux(:,uc_layidx:Nlay),2)./cntrlen,'LineWidth',linewidth);
hold on;
plot(y_avg/1000,rho0*sum(circ_meanMomFlux(:,uc_layidx:Nlay),2)./cntrlen,'LineWidth',linewidth);
plot(y_avg/1000,rho0*sum(circ_eddyMomFlux(:,uc_layidx:Nlay),2)./cntrlen,'LineWidth',linewidth);
hold off;
legend('Total','Mean','Eddy','Location','NorthWest');
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
title(['Along-isobath advective momentum forcing, layers ',num2str(uc_layidx),'-',num2str(Nlay)]);
if (write_figs)
  print('-dpng','-r150',fullfile(figdir,'IsobathAdvForcing_Undercurrent.png'));
end