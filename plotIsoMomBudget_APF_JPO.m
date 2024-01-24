%%%
%%% plotIsoMomBudget_APF_JPO.m
%%%
%%% Plots diagnostics from time-averaged momentum budget in isobath coordinates for our
%%% undercurrent run with an APF, for our JPO paper.
%%%

%%% Location of runs on file system
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';

%%% Select model configuration and parameters
grid_size = 256; %%% Default 128
wind_stress = 0.05; %%% Default 0.05
rand_force = 0; %%% Default 0.75
num_canyons = 4; %%% Default 4
amp_canyons = 25; %%% Default 25
max_slope = 0.15; %%% Default 0.15
sb_width = 5; %%% Default 5
baro_force = 0.025; %%% Default 0
drag_coeff = 2; %%% Default 2

%%% Plotting options
linewidth = 2;
axlim_zon = 7e-2;
axlim_iso = 5e-2;
ylim_zon = [0 350];
coord_name = 'bathy';
switch (coord_name)
  case 'bathy'
    d_levs = [110 150 200 400 1000 2000 2500 3000 3500];
  case 'psi'
    d_levs = [0.1 0.5 1 3 5];
end
fontsize = 14;
axpos = zeros(3,4);
axpos(1,:) = [0.145 0.69 .80 .26];
axpos(2,:) = [0.145 0.37 .80 .26];
axpos(3,:) = [0.145 0.05 .80 .26];
axlabels = {'(a)','(b)','(c)'};
lab_size = [0.05 0.03];
rho0 = 1000;
markersize = 40;
framepos = [382   100   520   1040];

figure(111);
clf;
set(gcf,'Position',framepos);
defaultcolor = get(gca,'Colororder');

%%% Load diagnostics for this run
config = 'wind';
uc_layidx = 3;
run_name = constructRunName (config,grid_size,wind_stress, ...
          rand_force,num_canyons,amp_canyons,max_slope,sb_width,baro_force,drag_coeff);
loadParams;
load(fullfile('products',[run_name,'_momBalance.mat']));

%%% For along-isobath plots
ylim_iso = [y_avg(1)/1000 100];
d_idx = zeros(1,length(d_levs));
d_labels = cell(1,length(d_levs));
for m=1:length(d_levs)
  d_idx(m) = find(abs(dd-d_levs(m))<1e-15);
  d_labels{m} = num2str(d_levs(m));
end

UC_yidx1 = find(dd>150,1,'First');
UC_yidx2 = find(dd>600,1,'First');


subplot('Position',axpos(1,:));
plot(y_avg/1000,sum(circ_windStress,2)./cntrlen,'LineWidth',linewidth,'Color',defaultcolororder(1,:));
hold on;
plot(y_avg/1000,sum(circ_Montgomery,2)./cntrlen,'LineWidth',2,'Color',defaultcolororder(2,:),'LineStyle',':');
plot(y_avg/1000,sum(circ_advection,2)./cntrlen,'LineWidth',linewidth,'Color',defaultcolororder(3,:));
plot(y_avg/1000,sum(circ_quadBotDrag,2)./cntrlen,'LineWidth',linewidth,'Color',defaultcolororder(4,:));
plot(y_avg/1000,sum(circ_hypervisc,2)./cntrlen,'LineWidth',linewidth,'Color',defaultcolororder(5,:));
plot(y_avg/1000,sum(circ_baroForcing,2)./cntrlen,'LineWidth',linewidth,'Color',defaultcolororder(6,:));
plot(y_avg/1000,sum(circ_total,2)./cntrlen,'LineWidth',linewidth,'Color',[.7 .7 .7]);
area([y_avg(UC_yidx1) y_avg(UC_yidx2)]/1000,[axlim_iso axlim_iso],-axlim_iso,'FaceColor','k','FaceAlpha',0.1);
hold off;
set(gca,'XDir','reverse');
ylabel('N/m^2');
set(gca,'YLim',[-axlim_iso axlim_iso]);
set(gca,'YTick',[-axlim_iso -axlim_iso/2 0 axlim_iso/2 axlim_iso]);
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
text(99,-axlim_iso*0.9,'Full depth','FontSize',fontsize+2);

subplot('Position',axpos(2,:));
plot(y_avg/1000,sum(circ_windStress(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth,'Color',defaultcolororder(1,:));
hold on;
plot(y_avg/1000,sum(circ_Montgomery(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth,'Color',defaultcolororder(2,:),'LineStyle','--');
plot(y_avg/1000,sum(circ_advection(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth,'Color',defaultcolororder(3,:));
plot(y_avg/1000,sum(circ_quadBotDrag(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth,'Color',defaultcolororder(4,:));
plot(y_avg/1000,sum(circ_hypervisc(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth,'Color',defaultcolororder(5,:));
plot(y_avg/1000,sum(circ_baroForcing(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth,'Color',defaultcolororder(6,:));
plot(y_avg/1000,sum(circ_total(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth,'Color',[.7 .7 .7]);
area([y_avg(UC_yidx1) y_avg(UC_yidx2)]/1000,[axlim_iso axlim_iso],-axlim_iso,'FaceColor','k','FaceAlpha',0.1);
hold off;
set(gca,'XDir','reverse');
ylabel('N/m^2');
set(gca,'YLim',[-axlim_iso axlim_iso]);
set(gca,'YTick',[-axlim_iso -axlim_iso/2 0 axlim_iso/2 axlim_iso]);
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
%   xlabel(ax2,'Isobath (m)');
box off;  
text(99,-axlim_iso*0.9,'Upper ocean','FontSize',fontsize+2);

subplot('Position',axpos(3,:));
plot(y_avg/1000,0*dd,'k--');
xlabel('Mean isobath offshore distance (km)');
set(gca,'XDir','reverse');
ylabel('N/m^2');
set(gca,'YLim',[-axlim_iso axlim_iso]);
set(gca,'YTick',[-axlim_iso -axlim_iso/2 0 axlim_iso/2 axlim_iso]);
set(gca,'XLim',ylim_iso);
set(gca,'FontSize',fontsize);
box off;
ax1 = gca;
ax2 = axes('Position',get(ax1,'Position'));
plot(y_avg/1000,sum(circ_windStress(:,uc_layidx:Nlay),2)./cntrlen,'LineWidth',linewidth,'Color',defaultcolororder(1,:));
hold on;
plot(y_avg/1000,-sum(circ_Montgomery(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth,'Color',defaultcolororder(2,:),'LineStyle','--');
plot(y_avg/1000,sum(circ_Montgomery(:,1:Nlay),2)./cntrlen,'LineWidth',linewidth,'Color',defaultcolororder(2,:),'LineStyle',':');
plot(y_avg/1000,sum(circ_advection(:,uc_layidx:Nlay),2)./cntrlen,'LineWidth',linewidth,'Color',defaultcolororder(3,:));
plot(y_avg/1000,sum(circ_quadBotDrag(:,uc_layidx:Nlay),2)./cntrlen,'LineWidth',linewidth,'Color',defaultcolororder(4,:));
plot(y_avg/1000,sum(circ_hypervisc(:,uc_layidx:Nlay),2)./cntrlen,'LineWidth',linewidth,'Color',defaultcolororder(5,:));
plot(y_avg/1000,sum(circ_baroForcing(:,uc_layidx:Nlay),2)./cntrlen,'LineWidth',linewidth,'Color',defaultcolororder(6,:));
plot(y_avg/1000,sum(circ_total(:,uc_layidx:Nlay),2)./cntrlen,'LineWidth',linewidth,'Color',[.7 .7 .7]);
area([y_avg(UC_yidx1) y_avg(UC_yidx2)]/1000,[axlim_iso axlim_iso],-axlim_iso,'FaceColor','k','FaceAlpha',0.1);
hold off;
handle = legend('Wind','Interfacial form stress','Topographic form stress','Advection','Friction','Viscosity','APF','Total','Location','SouthEast');
set(handle,'Position',[0.6250    0.2135    0.3452    0.1226]);
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
%   xlabel(ax2,'Isobath (m)');
box off;  
text(99,-axlim_iso*0.9,'Undercurrent','FontSize',fontsize+2);

%%% Add axis labels
for cntr = 1:size(axpos,1)
  annotation('textbox',[axpos(cntr,1)-0.12 axpos(cntr,2)-0.05 lab_size],'String',axlabels{cntr},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
end