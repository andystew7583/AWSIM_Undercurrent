%%%
%%% plotMomBudget.m
%%%
%%% Plots diagnostics from time-averaged momentum budget for our
%%% undercurrent runs.
%%%

%%% Location of runs on file system
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';

%%% Select model configuration and parameters
config = 'wind';
grid_size = 256; %%% Default 128
wind_stress = 0.05; %%% Default 0.05
rand_force = 0.75; %%% Default 0.75
num_canyons = 4; %%% Default 4
amp_canyons = 25; %%% Default 25
max_slope = 0.15; %%% Default 0.15
sb_width = 5; %%% Default 5
baro_force = 0.025; %%% Default 0
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
load(fullfile('products',[run_name,'_momBalance.mat']));

%%% For writing figures
write_figs = true;
figdir = fullfile('plots',run_name);
mkdir(figdir);

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