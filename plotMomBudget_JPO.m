%%%
%%% plotMomBudget_JPO.m
%%%
%%% Plots diagnostics from time-averaged momentum budget for our
%%% undercurrent runs, for our JPO paper.
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
drag_coeff = 3; %%% Default 2

%%% Plotting options
linewidth = 1.5;
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
axpos(1,:) = [0.075 0.7 .40 .26];
axpos(2,:) = [0.575 0.7 .40 .26];
axpos(3,:) = [0.075 0.38 .40 .26];
axpos(4,:) = [0.575 0.38 .40 .26];
axpos(5,:) = [0.075 0.06 .40 .26];
axpos(6,:) = [0.575 0.06 .40 .26];
axlabels = {'(a)','(b)','(c)','(d)','(e)','(f)'};
lab_size = [0.05 0.03];
rho0 = 1000;
markersize = 40;
framepos = [382   100   920   1040];
defaultcolororder = get(gca,'ColorOrder');

figure(106);
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
  load(fullfile('products',[run_name,'_momBalance.mat']));

  subplot('Position',axpos(1+i-1,:));
  plot(yy_h/1000,mean(sum(UMom_windStress,3),1),'LineWidth',linewidth,'Color',defaultcolororder(1,:));
  hold on;  
  plot(yy_h/1000,mean(sum(UMom_Montgomery,3),1),'LineWidth',linewidth,'Color',defaultcolororder(2,:),'LineStyle',':');
  plot(yy_h/1000,mean(sum(UMom_advection,3),1),'LineWidth',linewidth,'Color',defaultcolororder(3,:));
  plot(yy_h/1000,mean(sum(UMom_quadBotDrag,3),1),'LineWidth',linewidth,'Color',defaultcolororder(4,:));
  plot(yy_h/1000,mean(sum(UMom_hypervisc,3),1),'LineWidth',linewidth,'Color',defaultcolororder(5,:));
  plot(yy_h/1000,mean(sum(UMom_randomForcing,3),1),'LineWidth',linewidth,'Color',defaultcolororder(6,:));
  plot(yy_h/1000,mean(sum(UMom_total,3),1),'LineWidth',linewidth,'Color',[.7 .7 .7]);  
  hold off 
  set(gca,'XDir','reverse');  
%   if (i == 2)
%     handle = legend('Wind','Advection','Friction','Viscosity','Random forcing','Total','Location','NorthEast');
%     set(handle,'Position',[0.828804347826087 0.849519231322909 0.142934782608696 0.107692307692308]);
%   end
  if (i==1)
    ylabel('N/m^2');
  end
  set(gca,'YLim',[-axlim_zon axlim_zon]);
  set(gca,'XLim',ylim_zon);
  set(gca,'FontSize',fontsize);
  text(345,-axlim_zon*0.9,'Full depth','FontSize',fontsize+2);
  if (i == 1)
    title('Wind-forced','FontSize',fontsize+2);
  else
    title('Randomly forced','FontSize',fontsize+2);
  end

  subplot('Position',axpos(3+i-1,:));
  plot(yy_h/1000,mean(sum(UMom_windStress(:,:,1:uc_layidx-1),3),1),'LineWidth',linewidth,'Color',defaultcolororder(1,:));
  hold on;
  plot(yy_h/1000,mean(sum(UMom_Montgomery(:,:,1:uc_layidx-1),3),1),'LineWidth',linewidth,'LineStyle','--','Color',defaultcolororder(2,:));
  plot(yy_h/1000,mean(sum(UMom_advection(:,:,1:uc_layidx-1),3),1),'LineWidth',linewidth,'Color',defaultcolororder(3,:));
  plot(yy_h/1000,mean(sum(UMom_quadBotDrag(:,:,1:uc_layidx-1),3),1),'LineWidth',linewidth,'Color',defaultcolororder(4,:));
  plot(yy_h/1000,mean(sum(UMom_hypervisc(:,:,1:uc_layidx-1),3),1),'LineWidth',linewidth,'Color',defaultcolororder(5,:));
  plot(yy_h/1000,mean(sum(UMom_randomForcing(:,:,1:uc_layidx-1),3),1),'LineWidth',linewidth,'Color',defaultcolororder(6,:));
  plot(yy_h/1000,mean(sum(UMom_total(:,:,1:uc_layidx-1),3),1),'LineWidth',linewidth,'Color',[.7 .7 .7]);
  hold off  
  set(gca,'XDir','reverse');
  if (i==1)
    ylabel('N/m^2');
  end
  set(gca,'YLim',[-axlim_zon axlim_zon]);
  set(gca,'XLim',ylim_zon);
  set(gca,'FontSize',fontsize);  
  text(345,-axlim_zon*0.9,'Upper ocean','FontSize',fontsize+2);

  subplot('Position',axpos(5+i-1,:));
  plot(yy_h/1000,mean(sum(UMom_windStress(:,:,uc_layidx:Nlay),3),1),'LineWidth',linewidth,'Color',defaultcolororder(1,:));  
  hold on;    
  plot(yy_h/1000,mean(sum(-UMom_Montgomery(:,:,1:uc_layidx-1),3),1),'--','LineWidth',linewidth,'Color',defaultcolororder(2,:));
  plot(yy_h/1000,mean(sum(UMom_Montgomery(:,:,1:end),3),1),':','LineWidth',linewidth,'Color',defaultcolororder(2,:));
%   plot(yy_h/1000,mean(sum(UMom_Montgomery(:,:,uc_layidx:Nlay),3),1),'LineWidth',linewidth);
  plot(yy_h/1000,mean(sum(UMom_advection(:,:,uc_layidx:Nlay),3),1),'LineWidth',linewidth,'Color',defaultcolororder(3,:));
  plot(yy_h/1000,mean(sum(UMom_quadBotDrag(:,:,uc_layidx:Nlay),3),1),'LineWidth',linewidth,'Color',defaultcolororder(4,:));
  plot(yy_h/1000,mean(sum(UMom_hypervisc(:,:,uc_layidx:Nlay),3),1),'LineWidth',linewidth,'Color',defaultcolororder(5,:));
  plot(yy_h/1000,mean(sum(UMom_randomForcing(:,:,uc_layidx:Nlay),3),1),'LineWidth',linewidth,'Color',defaultcolororder(6,:));
  plot(yy_h/1000,mean(sum(UMom_total(:,:,uc_layidx:Nlay),3),1),'LineWidth',linewidth,'Color',[.7 .7 .7]);  
  hold off
  xlabel('Offshore distance (km)');  
  set(gca,'XDir','reverse');
  if (i == 2)
    handle = legend('Wind','Interfacial form stress','Topographic form stress','Advection','Friction','Viscosity','Random forcing','Total','Location','NorthEast');    
    set(handle,'Position',[0.791847826086958 0.530288461538465 0.195108695652174 0.122596153846154]);
  end
  if (i==1)
    ylabel('N/m^2');
  end
  set(gca,'YLim',[-axlim_zon axlim_zon]);
  set(gca,'XLim',ylim_zon);
  set(gca,'FontSize',fontsize); 
  text(345,-axlim_zon*0.9,'Undercurrent','FontSize',fontsize+2);

end


%%% Add axis labels
for cntr = 1:size(axpos,1)
  annotation('textbox',[axpos(cntr,1)-0.07 axpos(cntr,2)-0.05 lab_size],'String',axlabels{cntr},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
end