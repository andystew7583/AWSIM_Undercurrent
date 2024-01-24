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
linewidth = 1.5;
axlim_zon = 7e-2;
axlim_iso = 5e-2;
ylim_zon = [0 350];
d_levs = [110 150 200 400 1000 2000 2500 3000 3500];
d_idx = zeros(1,length(d_levs));
d_labels = cell(1,length(d_levs));
fontsize = 14;
axpos = zeros(4,4);
axpos(1,:) = [0.075 0.73 .40 .23];
axpos(2,:) = [0.535 0.73 .40 .23];
axpos(3,:) = [0.075 0.45 .40 .23];
axpos(4,:) = [0.535 0.45 .40 .23];
axpos(5,:) = [0.075 0.06 .40 .33];
axpos(6,:) = [0.535 0.06 .40 .33];
cbpos = zeros(3,4);
cbpos(1,:) = [0.95 0.73 0.01 0.23];
cbpos(2,:) = [0.95 0.45 0.01 0.23];
cbpos(3,:) = [0.95 0.06 0.01 0.33];
axlabels = {'(a)','(b)','(c)','(d)','(e)','(f)'};
lab_size = [0.05 0.03];
rho0 = 1000;
markersize = 40;
framepos = [382   100   920   850];

figure(104);
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
  ylim_iso = [y_avg(1)/1000 100];
  for m=1:length(d_levs)
    d_idx(m) = find(abs(dd-d_levs(m))<1e-15);
    d_labels{m} = num2str(d_levs(m));
  end
  h_w = 0.5*(h(1:Nx,:,:)+h([Nx 1:Nx-1],:,:));   
  u_twa = hu./h_w;  
  
  %%% Create grid for slice contour plots
  eps = 1; %%% Small distance by which to perturb plotting points away from actual isopycnal elevations
  slice_idx = Nx/4;
  % slice_idx = 3*Nx/8;
  ZZ_slice = zeros(Ny,2*Nlay);
  UU_slice = zeros(Ny,2*Nlay);
  for k=1:Nlay
    ZZ_slice(:,2*k-1) = squeeze(e(slice_idx,:,k))-eps;
    ZZ_slice(:,2*k) = squeeze(e(slice_idx,:,k+1))+eps;
    UU_slice(:,2*k-1) = squeeze(u_twa(slice_idx,:,k));
    UU_slice(:,2*k) = squeeze(u_twa(slice_idx,:,k));
  end
  YY_slice = repmat(yy_h',[1 2*Nlay]);
  
  subplot('Position',axpos(1+i-1,:));
  pcolor(XX_u/1000,YY_u/1000,u_twa(:,:,1))
  shading interp
  hold on;
  [C,handle] = contour(XX_h/1000,YY_h/1000,-hhb,[110 300 1500 3000 3500 4000],'EdgeColor','k');
  clabel(C,handle);
  plot([xx_u(slice_idx) xx_u(slice_idx)]/1000,[0 yy_u(end)/1000],'k--','LineWidth',3);
  hold off  
  colormap(gca,cmocean('balance',40))
  caxis([-.4 .4]);
  set(gca,'FontSize',fontsize);
  if (i == 1)
    ylabel('Offshore distance (km)');
  end
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

  subplot('Position',axpos(3+i-1,:));
  pcolor(XX_u/1000,YY_u/1000,u_twa(:,:,uc_layidx))
  shading interp
  hold on;
  [C,handle] = contour(XX_h/1000,YY_h/1000,-hhb,[110 300 1500 3000 3500 4000],'EdgeColor','k');
  clabel(C,handle);
  plot([xx_u(slice_idx) xx_u(slice_idx)]/1000,[0 yy_u(end)/1000],'k--','LineWidth',3);
  hold off  
  colormap(gca,cmocean('balance',20))
  caxis([-.2 .2]);
  set(gca,'FontSize',fontsize);
  xlabel('Alongshore distance (km)');
  if (i == 1)
    ylabel('Offshore distance (km)');
  end
  if (i == 2)
    handle = colorbar;
    set(handle,'Position',cbpos(2,:));
    set(get(handle,'Title'),'String','m/s')
  end

  subplot('Position',axpos(5+i-1,:));
  pcolor(YY_slice/1000,-ZZ_slice,UU_slice);
  shading interp
  hold on;
  plot(yy_h/1000,-squeeze(e(slice_idx,:,:)),'Color',[.7 .7 .7]);  
  hold off;
  colormap(gca,cmocean('balance',30))
  caxis([-0.3 0.3]);
  xlabel('Offshore distance (km)');
  if (i == 1)    
    ylabel('Depth (m)');
  end
  set(gca,'XDir','reverse');
  set(gca,'YDir','reverse');
  set(gca,'XLim',ylim_zon);
  set(gca,'FontSize',fontsize);
  if (i == 2)
    handle = colorbar;
    set(handle,'Position',cbpos(3,:));
    set(get(handle,'Title'),'String','m/s')
  end
  set(gca,'Color',[.8 .8 .8]);



end




%%% Add axis labels
for cntr = 1:size(axpos,1)
  annotation('textbox',[axpos(cntr,1)-0.07 axpos(cntr,2)-0.05 lab_size],'String',axlabels{cntr},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
end