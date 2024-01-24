%%%
%%% plotAlongIsobathFlow_JPO.m
%%%
%%% Plots plots along-isobath-averaged flow diagnostics for our JPO paper.
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
axpos = zeros(4,4);
axpos(1,:) = [0.075 0.55 .39 .41];
axpos(2,:) = [0.535 0.55 .39 .41];
axpos(3,:) = [0.075 0.07 .39 .41];
axpos(4,:) = [0.535 0.07 .39 .41];
cbpos = zeros(2,4);
cbpos(1,:) = [0.95 0.55 0.01 0.41];
cbpos(2,:) = [0.95 0.07 0.01 0.41];
axlabels = {'(a)','(b)','(c)','(d)'};
lab_size = [0.05 0.03];
rho0 = 1000;
markersize = 40;
framepos = [382   100   920   720];

figure(105);
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

  %%% Create grid for along-isobath averaged contour plots
  eps = 1; %%% Small distance by which to perturb plotting points away from actual isopycnal elevations
  ZZ_avg = zeros(Nd,2*Nlay);
  UU_avg = zeros(Nd,2*Nlay);
  EE_avg = zeros(Nd,2*Nlay);
  for k=1:Nlay
    ZZ_avg(:,2*k-1) = e_avg(:,k)-eps;
    ZZ_avg(:,2*k) = e_avg(:,k+1)+eps;
    UU_avg(:,2*k-1) = circ_twa(:,k);
    UU_avg(:,2*k) = circ_twa(:,k);
    EE_avg(:,2*k-1) = EKE(:,k); %%% TODO remove once EKE calculation is updated
    EE_avg(:,2*k) = EKE(:,k);
  end
  YY_avg = repmat(y_avg,[1 2*Nlay]);
  
  

  subplot('Position',axpos(1+i-1,:));
  pcolor(YY_avg/1000,-ZZ_avg,UU_avg);
  shading interp
  hold on;
  plot(y_avg/1000,-e_avg,'Color',[.7 .7 .7]);
  plot(
  hold off;
  colormap(gca,cmocean('balance',40)); 
  caxis([-.1 .1])
  if (i == 1)
    ylabel('Depth (m)');
  end
  set(gca,'XLim',[y_avg(1)/1000 100]);
  set(gca,'YLim',[0 3000]);
  set(gca,'XDir','reverse');
  set(gca,'YDir','reverse');
  set(gca,'FontSize',fontsize);
  if (i == 2)
    handle = colorbar;
    set(handle,'Position',cbpos(1,:));
    set(get(handle,'Title'),'String','m/s')
  end
  if (i == 1)
    title('Wind-forced','FontSize',fontsize+2);
  else
    title('Randomly forced','FontSize',fontsize+2);
  end
  set(gca,'Color',[.8 .8 .8]); 
  text(78,2800,'Along-slope flow','FontSize',fontsize+2);
  
  subplot('Position',axpos(3+i-1,:));
  pcolor(YY_avg/1000,-ZZ_avg,log10(EE_avg));
  shading interp
  hold on;
  plot(y_avg/1000,-e_avg,'Color',[.7 .7 .7]);
  hold off;
  colormap(gca,cmocean('amp',25));  
  caxis([-4 -1.5])
  xlabel('Mean isobath offshore distance (km)');
  if (i == 1)
    ylabel('Depth (m)');
  end
  set(gca,'XLim',[y_avg(1)/1000 100]);
  set(gca,'YLim',[0 3000]);  
  set(gca,'XDir','reverse');
  set(gca,'YDir','reverse');
  set(gca,'FontSize',fontsize);  
  set(gca,'Color',[.8 .8 .8]);  
  text(83,2800,'Eddy kinetic energy','FontSize',fontsize+2);
  if (i == 2)
    handle = colorbar;
    set(handle,'Position',cbpos(2,:));
    set(get(handle,'Title'),'String','m^2/s^2')
    set(handle,'YTick',[-4:1:-2]);
    set(handle,'YTickLabel',{'10^-^4','10^-^3','10^-^2'});
  end
  
end

%%% Add panel labels
for cntr = 1:size(axpos,1)
  annotation('textbox',[axpos(cntr,1)-0.07 axpos(cntr,2)-0.05 lab_size],'String',axlabels{cntr},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
end