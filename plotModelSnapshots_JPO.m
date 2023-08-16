%%%
%%% plotModelSnapshots_JPO.m
%%%
%%% Plots 3D renderings of instantaneous model snapshots for our JPO paper.
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
axpos = zeros(2,4);
axpos(1,:) = [0.1 0.18 .35 .78];
axpos(2,:) = [0.6 0.18 .35 .78];
cbpos = [0.4 0.04 0.2 0.01];
axlabels = {'(a)','(b)'};
lab_size = [0.05 0.03];
rho0 = 1000;
markersize = 40;
framepos = [382   100   920   630];
n = 20; %%% Select output iteration number

figure(103);
clf;
set(gcf,'Position',framepos);

configs = {'wind','rand'};
for i=1:2
  
  %%% Load diagnostics for this run
  config = configs{i};  
  run_name = constructRunName (config,grid_size,wind_stress, ...
            rand_force,num_canyons,amp_canyons,max_slope,sb_width,baro_force,drag_coeff);
  loadParams;

  %%% Calculate layer surface height
  eta = zeros(Nx,Ny,Nlay+1);
  eta(:,:,Nlay+1) = hhb;
  hh = zeros(Nx,Ny,Nlay);
  zeta = zeros(Nx+1,Ny+1,Nlay);            
  for k=Nlay:-1:1

    %%% Load kth layer thickness
    data_file = fullfile(dirpath,[OUTN_H,num2str(k-1),'_n=',num2str(n),'.dat']);
    hh(:,:,k) = readOutputFile(data_file,Nx,Ny);    

    %%% Add layer thickness to eta
    eta(:,:,k) = eta(:,:,k+1) + hh(:,:,k);

    %%% Calculate the relative vorticity
    data_file = fullfile(dirpath,[OUTN_U,num2str(k-1),'_n=',num2str(n),'.dat']);
    uu = readOutputFile(data_file,Nx,Ny);
    data_file = fullfile(dirpath,[OUTN_V,num2str(k-1),'_n=',num2str(n),'.dat']);
    vv = readOutputFile(data_file,Nx,Ny);           
    zeta(2:Nx,2:Ny,k) = (uu(1:Nx-1,1:Ny-1)-uu(1:Nx-1,2:Ny)) / dy + (vv(2:Nx,2:Ny)-vv(1:Nx-1,2:Ny)) / dx;                      

  end

  %%% Correct boundaries
  if (~useWallNS)
    zeta(:,Ny+1,:) = zeta(:,1,:);
  end
  zeta(Nx+1,:,:) = zeta(1,:,:);        


  %%% Make the plot
  subplot('Position',axpos(1+i-1,:));
  pbaspect([2,1,1]);
  eta_plot = eta(:,:,2);
  Ro = zeta(:,:,1)./(2*Omega_z);
  p = surface(XX_h/1000,YY_h/1000,eta_plot,Ro);
%   p = surface(XX_h/1000,YY_h/1000,eta_plot,sqrt(uu(:,:,1).^2+vv(:,:,1).^2));
  p.FaceColor = 'texturemap';
  colormap(cmocean('balance'));
  if (i == 1)
    handle = colorbar;
    set(handle,'Position',cbpos);
    set(handle,'Orientation','Horizontal');
    set(get(handle,'Title'),'String','Vortex Rossby number');
  end  
  caxis([-.5 .5]);
  p.EdgeColor = 'none';         
  alpha(p,0.7);
  hold on;
  for k = 3:Nlay
    eta_plot = eta(:,:,k);
    p = surface(XX_h/1000,YY_h/1000,eta_plot);
    p.FaceColor = [48 129 238]/256;
    p.EdgeColor = 'none';        
    alpha(p,0.5);          
  end
  p = surface(XX_h/1000,YY_h/1000,eta(:,:,Nlay+1));
  p.FaceColor = [139,69,19]/256;
  p.EdgeColor = 'none';          
  hold off;        
  view(100,30);
  lighting gouraud;
  camlight(90,60);
%         camlight(-45,60);
%         camlight(180,60);        
  if (i == 1)
    title('Wind-forced','FontSize',fontsize+2);
  else
    title('Randomly forced','FontSize',fontsize+2);
  end
  xlabel('x (km)');
  ylabel('y (km)');
  zlabel('z (m)');
  set(gca,'FontSize',16);
  set(gcf,'Color','w');        
  set(gca,'ZLim',[-4650 0]);

end


%%% Add axis labels
for cntr = 1:size(axpos,1)
  annotation('textbox',[axpos(cntr,1)-0.07 axpos(cntr,2)-0.05 lab_size],'String',axlabels{cntr},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
end
