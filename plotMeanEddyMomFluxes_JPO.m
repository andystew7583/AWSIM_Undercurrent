%%%
%%% plotMeanFlow.m
%%%
%%% Plots diagnostics from time-averaged mean flow for our
%%% undercurrent runs.
%%%

%%% Location of runs on file system
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';

%%% Select model configuration and parameters
grid_size = 256; %%% Default 128
wind_stress = 0.05; %%% Default 0.05
rand_force = 0.75; %%% Default 0.75
num_canyons = 4; %%% Default 4
amp_canyons = 25; %%% Default 25
max_slope = 0.15; %%% Default 0.15
sb_width = 5; %%% Default 5
baro_force = 0; %%% Default 0
drag_coeff = 2; %%% Default
       
%%% Plotting options
linewidth = 2;
axlim_zon = 7e-2;
axlim_iso = 5e-2;
ylim_iso = [y_avg(1)/1000 100];
ylim_zon = [0 350];
switch (coord_name)
  case 'bathy'
    d_levs = [110 150 200 400 1000 2000 2500 3000 3500];
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
axpos = zeros(6,4);
axpos(1,:) = [0.075 0.68 .40 .25];
axpos(2,:) = [0.575 0.68 .40 .25];
axpos(3,:) = [0.075 0.37 .40 .25];
axpos(4,:) = [0.575 0.37 .40 .25];
axpos(5,:) = [0.075 0.06 .40 .25];
axpos(6,:) = [0.575 0.06 .40 .25];
axlabels = {'(a)','(b)','(c)','(d)','(e)','(f)'};
lab_size = [0.05 0.03];
rho0 = 1000;
markersize = 40;
framepos = [382   100   920   1040];

figure(108);
clf;
set(gcf,'Position',framepos);

configs = {'wind','rand'};
uc_layidxs = [3 2];
for i=1:2
  
  %%% Load diagnostics for this run
  config = configs{i};
  uc_layidx = uc_layidxs(i);
  run_name = constructRunName (config,grid_size,wind_stress, ...
            rand_force,num_canyons,amp_canyons,max_slope,sb_width,baro_force,drag_coeff);
  loadParams;
  load(fullfile('products',[run_name,'_meanFlow.mat']));
  
  UC_yidx1 = find(dd>150,1,'First');
  UC_yidx2 = find(dd>600,1,'First');


  %%% For momentum balance plots
  axlim_iso = 5e-2;

  subplot('Position',axpos(1+i-1,:));
  plot(y_avg/1000,rho0*sum(circ_totalMomFlux,2)./cntrlen,'LineWidth',linewidth);
  hold on;
  plot(y_avg/1000,rho0*sum(circ_meanMomFlux,2)./cntrlen,'LineWidth',linewidth);
  plot(y_avg/1000,rho0*sum(circ_eddyMomFlux,2)./cntrlen,'LineWidth',linewidth);
  area([y_avg(UC_yidx1) y_avg(UC_yidx2)]/1000,[axlim_iso axlim_iso],-axlim_iso,'FaceColor','k','FaceAlpha',0.1);
  hold off; 
  set(gca,'XDir','reverse');
  if (i == 1)
    ylabel('N/m^2');
  end
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
  if (i == 1)
    title('Wind-forced','FontSize',fontsize+2);
  else
    title('Randomly forced','FontSize',fontsize+2);
  end
  text(99,-axlim_iso*0.9,'Full depth','FontSize',fontsize+2);

  subplot('Position',axpos(3+i-1,:));
  plot(y_avg/1000,rho0*sum(circ_totalMomFlux(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth);
  hold on;
  plot(y_avg/1000,rho0*sum(circ_meanMomFlux(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth);
  plot(y_avg/1000,rho0*sum(circ_eddyMomFlux(:,1:uc_layidx-1),2)./cntrlen,'LineWidth',linewidth);
  area([y_avg(UC_yidx1) y_avg(UC_yidx2)]/1000,[axlim_iso axlim_iso],-axlim_iso,'FaceColor','k','FaceAlpha',0.1);
  hold off;
  if (i == 2)
    legend('Total','Mean','Eddy','Location','NorthWest'); 
  end
  set(gca,'XDir','reverse');
  if (i == 1)
    ylabel('N/m^2');
  end
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
  box off; 
  text(99,-axlim_iso*0.9,'Upper ocean','FontSize',fontsize+2);

  subplot('Position',axpos(5+i-1,:));
  plot(y_avg/1000,rho0*sum(circ_totalMomFlux(:,uc_layidx:Nlay),2)./cntrlen,'LineWidth',linewidth);
  hold on;
  plot(y_avg/1000,rho0*sum(circ_meanMomFlux(:,uc_layidx:Nlay),2)./cntrlen,'LineWidth',linewidth);
  plot(y_avg/1000,rho0*sum(circ_eddyMomFlux(:,uc_layidx:Nlay),2)./cntrlen,'LineWidth',linewidth);
  area([y_avg(UC_yidx1) y_avg(UC_yidx2)]/1000,[axlim_iso axlim_iso],-axlim_iso,'FaceColor','k','FaceAlpha',0.1);
  hold off;  
  xlabel('Mean isobath offshore distance (km)');
  set(gca,'XDir','reverse');
  if (i == 1)
    ylabel('N/m^2');
  end
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
  box off; 
  text(99,-axlim_iso*0.9,'Undercurrent','FontSize',fontsize+2);
  
end

%%% Add panel labels
for cntr = 1:size(axpos,1)
  annotation('textbox',[axpos(cntr,1)-0.07 axpos(cntr,2)-0.05 lab_size],'String',axlabels{cntr},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
end