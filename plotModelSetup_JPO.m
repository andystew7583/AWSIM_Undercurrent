%%%
%%% plotModelSetup_JPO.m
%%%
%%% Plots wind forcing, instantaneous model state and bathymetry for our JPO paper.
%%%

%%% Location of runs on file system
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';

%%% Run to load
% run_name = 'undercurrent_LWF_tau0.05_Nc4_Yc25_Ny256';
run_name = 'undercurrent_LWF_tau0.05_Nc4_Yc25_Fbaro0.025_Ny256';

%%% Load parameters   
loadParams;
dirpath = fullfile(local_home_dir,run_name);
gtild = reshape(cumsum(gg),[Nlay 1 1]);

%%% Load constant parameters 
constants;   

%%% Physical parameters
f0 = 8e-5;                   %%% Coriolis parameter 
beta = 0e-11;               %%% Coriolis parameter gradient
Ly = 400*m1km;                %%% Domain length.   
Lx = 2*Ly;
H = 4000;                     %%% Ocean depth  
Cd = 2e-3;                %%% Quadratic bottom drag
tau0 = 0.05;                  %%% Wind stress maximum      
tRelax_u = 1*t1day;         %%% Sponge relaxation time scales
tRelax_v = 1*t1day;
tRelax_h = 1*t1day;          
Lrelax = 50*m1km;           %%% Sponge relaxation width

%%% Exponentials and linear slope with matched gradients
smax = max_slope; %%% Max slope steepness
Ws_shelf = sb_width*m1km; %%% Reference slope exponential width for continental shelf
Ws_deep = 20*m1km; %%% Reference slope exponential width for deep ocean
Ys_shelf = 75*m1km; %%% Latitude of center of continental slope
H_shelf = 100; %%% Continental shelf water column thickness  
DY_canyons = amp_canyons*m1km; %%% Meridional amplitude of slope canyons
N_canyons = num_canyons; %%% Number of slope canyons
X_canyons = Lx/4-Lx/N_canyons/2; %%% Longitude of first slope canyon
H_ridge = 0; %%% Meridional ridge height
W_ridge = 50*m1km; %%% Zonal width of meridional ridge
etab_lambda = 100*m1km;
etab_K_0 = 2*pi/etab_lambda; %%% Most energetic wavenumber
etab_W = etab_K_0/8; %%% Exponential width of energy band in wavenumber space

%%% Matrices defining spatially-dependent geometric parameters
Zbot = -H * ones(Nx,Ny); %%% Flat reference ocean floor depth    
Zshelf = -H_shelf * ones(Nx,Ny); %%% Continental shelf elevation
Ws_upper = Ws_shelf * ones(Nx,Ny); %%% Exponential width of upper slope
Ws_lower = Ws_deep * ones(Nx,Ny); %%% Exponential width of lower slope, reverts to Ws_shelf where ridge elevation is high
Delta_upper = Ws_upper * smax; %%% Exponential elevation change of upper slope
Delta_lower = Ws_lower * smax; %%% Exponential elevation change of lower slope
Zs_upper = (Zshelf - Delta_upper); %%% Elevation of upper/mid slope join
Zs_lower = (Zbot + Delta_lower); %%% Elevation of lower/mid slope join  
Ys_upper = Ys_shelf - DY_canyons*cos(2*pi*N_canyons*(XX_h-X_canyons)/Lx); %%% Latitude of upper/mid slope join
Ys_lower = Ys_upper + ((Zs_upper-Zs_lower)/smax); %%% Latitude of lower/mid slope join

%%% Construct the bathymetry by joining exponential upper and lower
%%% slopes with a linear mid-slope
etab = zeros(Nx,Ny);
idx_upper = YY_h < Ys_upper;
idx_mid = (YY_h >= Ys_upper) & (YY_h < Ys_lower);
idx_lower = YY_h > Ys_lower;    
etab_upper = Zshelf - Delta_upper .* exp((YY_h-Ys_upper)./Ws_upper);
etab_mid = Zs_upper - smax.*(YY_h-Ys_upper);
etab_lower = Zbot + Delta_lower .* exp(-(YY_h-Ys_lower)./Ws_lower);
etab(idx_upper) = etab_upper(idx_upper);
etab(idx_mid) = etab_mid(idx_mid);
etab(idx_lower) = etab_lower(idx_lower);

%%% Stratification parameters
rho0 = 1000;                  %%% Reference density 
dz = 1;                       %%% Discretization step for continuous density profile
zz = 0:-dz:-H;                %%% Grid for continuous density profile
hs = 400;                     %%% Exponential decay scale for density profile
hl = 40000;                   %%% Linear decay scale for density profile 
rhoBot = 1028;                %%% Sea floor density in continuous density profile
rhoSurf = 1024;               %%% Sea surface density in continuous density profile
rhoRange = rhoBot-rhoSurf;

%%% Remove redundant time dimension in wind stress vector
taux = squeeze(mean(taux,1));
tauy = squeeze(mean(tauy,1));

%%% Exponential/inear continuous density profile
rhoRef = rhoBot - rhoRange*(exp(zz/hs)+zz/hl-exp(-H/hs)+H/hl)/(1-exp(-H/hs)+H/hl); 

%%% Piecewise density profile for wind-forced cases
Hupper = 150; 
Hinner = 150;   
E0_wind = 0;
if (Nlay >= 1)
  E0_wind = [E0_wind -Hupper];
end
if (Nlay >= 2)
  E0_wind = [E0_wind -Hupper-Hinner*(1:Nlay-2)];
end  
E0_wind = [E0_wind (-H)];
H0 = E0_wind(1:Nlay)-E0_wind(2:Nlay+1);
rhoMean_wind = zeros(1,Nlay);
rhoProf_wind = 0*zz;
for k=1:Nlay
  idx = find((zz<=E0_wind(k))  & (zz>=E0_wind(k+1)));
  rhoMean_wind(k) = mean(rhoRef(idx));
  rhoProf_wind(idx) = rhoMean_wind(k);
end

%%% Piecewise density profile for randomly forced cases
Nlay = 4;
Hupper = 50;
Hinner = 200; 
E0_rand = 0;
if (Nlay >= 1)
  E0_rand = [E0_rand -Hupper];
end
if (Nlay >= 2)
  E0_rand = [E0_rand -Hupper-Hinner*(1:Nlay-2)];
end  
E0_rand = [E0_rand (-H)];
H0 = E0_rand(1:Nlay)-E0_rand(2:Nlay+1);
rhoMean_rand = zeros(1,Nlay);
rhoProf_rand = 0*zz;
for k=1:Nlay
  idx = find((zz<=E0_rand(k))  & (zz>=E0_rand(k+1)));
  rhoMean_rand(k) = mean(rhoRef(idx));
  rhoProf_rand(idx) = rhoMean_rand(k);
end 

%%% Random forcing parameters
useBaroclinicForcing = true; %%% Apply random forcing to near-surface layers only  
RF_tau = 864000; %%% Autocorrelation timescale  
lambdaK = 80*m1km; %%% Peak energy input lengthscale for random forcing
Ymin_RF = 0; %%% No random forcing below this latitude
Ymax_RF = Ly-Lrelax; %%% No random forcing above this latitude
Lramp_RF = 30*m1km; %%% Latitudinal width of region over which to ramp up the random forcing

  
%%% Generate forcing mask in spectral space
k = [0:1:Nx/2-1,-Nx/2:1:-1]; %%% Zonal wavenumber
K_xk = 2*pi.*(k)./Lx;
l = [0:1:Ny/2-1,-Ny/2:1:-1]; %%% Meridional wavenumber
K_yl = 2*pi.*(l)./Ly;
[K_ykl,K_xkl] = meshgrid(K_yl, K_xk);
K = sqrt(K_xkl.^2 + K_ykl.^2); %%% Absolute wavenumber 
K_0 = 2*pi/lambdaK; %%% Most strongly forced wavenumber
W = K_0/8; %%% Exponential width of energy band in wavenumber space   
RF_mask_fft = exp(-((K-K_0)/W).^2);
RF_mask_norm = sqrt(0.5*sum(RF_mask_fft(:)));

RF_F0_rot = 0.75;
RF_kk_recip = 1.0 ./ K;
RF_kk_recip(1,1) = 0;
RFrot_fft = complex(zeros(Nx,Ny));
RFamp_rot_fft = RF_F0_rot .* RF_kk_recip .* RF_mask_fft ./ RF_mask_norm;

%%% Evolve over many forcing decorrelation time scales using the actual
%%% time step for this run to get a representative forcing field
%%% N.B. This code is basically taken directly
deltaT = 65;
Tratio = deltaT/RF_tau;
for n=1:20000
  RF_phase = 2*pi*rand(Nx,Ny);
  RFrot_fft = RFrot_fft - Tratio.*RFrot_fft + sqrt(2*Tratio).*exp(1i*RF_phase).*RFamp_rot_fft;
end
RFrot = Nx*Ny*real(ifft2(RFrot_fft));

%%% Test that RMS forcing magnitude matches specified forcing amplitude
dRFrot_dx = (RFrot(1:Nx,1:Ny)-RFrot([Nx 1:Nx-1],1:Ny))/dx;
dRFrot_dy = (RFrot(1:Nx,1:Ny)-RFrot(1:Nx,[Ny 1:Ny-1]))/dy;
gradsq_dRFrot = dRFrot_dx.^2+dRFrot_dy.^2;
sqrt(mean(gradsq_dRFrot(:)))

%%% Generate forcing mask in real space
idx_south_zero = (yy_q(1:Ny) <= Ymin_RF);
idx_south_cos = (yy_q(1:Ny) <= Ymin_RF+Lramp_RF) & (yy_q(1:Ny) >= Ymin_RF);
idx_north_cos = (yy_q(1:Ny) >= Ymax_RF-Lramp_RF) & (yy_q(1:Ny) <= Ymax_RF);
idx_north_zero = (yy_q(1:Ny) >= Ymax_RF);
RF_mask_q = ones(Nlay,Nx,Ny);
if (~isempty(idx_south_zero))
  RF_mask_q(:,:,idx_south_zero) = 0;
end
if (~isempty(idx_south_cos))    
  RF_mask_q(:,:,idx_south_cos) = RF_mask_q(:,:,idx_south_cos) .* repmat(reshape(0.5.*(1+cos(pi*(YY_q(1:Nx,idx_south_cos)-Ymin_RF-Lramp_RF)/Lramp_RF)),[1 Nx sum(idx_south_cos)]),[Nlay 1 1]);
end
if (~isempty(idx_north_cos))
  RF_mask_q(:,:,idx_north_cos) = RF_mask_q(:,:,idx_north_cos) .* repmat(reshape(0.5.*(1+cos(pi*(YY_q(1:Nx,idx_north_cos)-Ymax_RF-Lramp_RF)/Lramp_RF)),[1 Nx sum(idx_north_cos)]),[Nlay 1 1]);
end
if (~isempty(idx_north_zero))
  RF_mask_q(:,:,idx_north_zero) = 0;
end 
RF_mask_q = RF_mask_q / H;
RF_mask_q = squeeze(RF_mask_q(1,:,:));

  
  
  
  
  
  

%%% Plotting options
scrsz = get(0,'ScreenSize');
framepos = [209   632   900   900];
fontsize = 14;
axpos = zeros(7,4);
axpos(1,:) = [0.06 0.71 0.25 0.26];
axpos(2,:) = [0.3925 0.71 0.53 0.26];
axpos(3,:) = [0.06 0.385 0.25 0.26];
axpos(4,:) = [0.3925 0.385 0.53 0.26];
axpos(5,:) = [0.06 0.06 0.25 0.26];
axpos(6,:) = [0.3925 0.06 0.25 0.26];
axpos(7,:) = [0.705 0.06 0.25 0.26];
cb_pos1 = [0.94 0.71 0.015 0.26];
cb_pos2 = [0.94 0.3925 0.015 0.26];
lab_size = [0.05 0.03];
axlabels = {'(a)','(b)','(c)','(d)','(e)','(f)','(g)'};


%%% Set up the frame
figure(102);
clf; 
set(gcf,'Position',framepos);
set(gcf,'Color','w');
cntr = 0;



%%% Plot wind stress
cntr = cntr + 1;
ax = axes('position',axpos(cntr,:));  
plot(rho0*taux(1,:),yy_h/m1km,'LineWidth',1.5);
hold on;
plot(rho0*mean(Fbaro_x.*-hhb,1),yy_h/m1km,'--','LineWidth',1.5);
plot([0 0],[0 Ly/m1km],'k:','LineWidth',0.5);
hold off;
axis([-0.06 0.06 0 Ly/m1km]);
ylabel('Offshore distance (km)');
xlabel('Alongshore stress (N/m^2)');
% title('Wind stress');
set(gca,'FontSize',fontsize);
legend('Wind stress','APF','Location','Northeast');


%%% Plot PV snapshot
cntr = cntr + 1;
ax = axes('position',axpos(cntr,:));  
pcolor(XX_h/1000,YY_h/1000,-hhb);
shading interp;
hold on;
[C,h]=contour(XX_h/1000,YY_h/1000,-hhb,[110 300 1500 3000 3500 4000],'EdgeColor',[.3 .3 .3],'LineWidth',1);
fill([0 Lx Lx 0 0]/m1km,[Ly-Lrelax Ly-Lrelax Ly Ly Ly-Lrelax]/m1km,ones(1,5),'FaceColor',[.5 .5 .5],'EdgeColor','None','FaceAlpha',0.5);
% plot([0 Lx Lx 0 0]/m1km,[Ly-Ln Ly-Ln Ly Ly Ly-Ln]/m1km,'k-');
clabel(C,h);
% plot([-Lx/2/m1km Lx/2/m1km Lx/2/m1km -Lx/2/m1km -Lx/2/m1km],[-Ly/2/m1km -Ly/2/m1km Ly/2/m1km Ly/2/m1km -Ly/2/m1km],'k-');
hold off;
% axis([-Lx/2/m1km Lx/2/m1km -Ly/2/m1km Ly/2/m1km]);
set(gca,'FontSize',fontsize);
% colormap(pmkmp(100,'Swtth'));
colormap(gca,flip(haxby));
cbhandle = colorbar;
set(cbhandle,'Position',cb_pos1);
set(get(cbhandle,'Title'),'String','m');
% caxis([min(pv(:)) max(pv(:))]);
% title('Model bathymetry');
text(Lx/2/m1km,(Ly-Lrelax/2)/m1km,'Sponge','FontSize',10);
title('Bathymetry');

%%% Plot random forcing shape in spectral space
cntr = cntr + 1;
ax = axes('position',axpos(cntr,:));  

KK = 1e-6:1e-6:max(K_xk);
K_0 = 2*pi/lambdaK; %%% Most strongly forced wavenumber
W = K_0/8; %%% Exponential width of energy band in wavenumber space   
semilogx(2*pi./KK/1000,exp(-((KK-K_0)/W).^2),'LineWidth',linewidth);
hold on;
semilogx(2*pi./KK/1000,exp(-((KK-etab_K_0)/etab_W).^2),'LineWidth',linewidth,'Color',[139,69,19]/256);
hold off;
axis([10 1000 0 1.3]);
ylabel('Normalized spectral amplitude');
xlabel('Wavelength (km)');
set(gca,'FontSize',fontsize);
legend('Random forcing','Bathymetric perturbations');

%%% Plot random forcing snapshot
cntr = cntr + 1;
ax = axes('position',axpos(cntr,:));  
pcolor(XX_h/1000,YY_h/1000,RFrot.*RF_mask_q);
shading interp;
hold on;
[C,h]=contour(XX_h/1000,YY_h/1000,-hhb,[110 300 1500 3000 3500 4000],'EdgeColor',[.3 .3 .3],'LineWidth',1);
fill([0 Lx Lx 0 0]/m1km,[Ly-Lrelax Ly-Lrelax Ly Ly Ly-Lrelax]/m1km,ones(1,5),'FaceColor',[.5 .5 .5],'EdgeColor','None','FaceAlpha',0.5);
clabel(C,h);
hold off;
set(gca,'FontSize',fontsize);
colormap(gca,cmocean('balance'));
cbhandle = colorbar;
set(cbhandle,'Position',cb_pos2);
set(get(cbhandle,'Title'),'String','N/m^2');
caxis([-10 10]);
xlabel('Alongshore distance (km)');
ylabel('Offshore distance (km)');
text(Lx/2/m1km,(Ly-Lrelax/2)/m1km,'Sponge','FontSize',10);
title('Random force function');

cntr = cntr + 1;
ax = axes('position',axpos(cntr,:));  
plot(yy_h/1000,etab(Nx/2,:)/1000,'-','LineWidth',linewidth,'Color',[139,69,19]/256);
hold on;
plot(yy_h/1000,etab(3*Nx/8,:)/1000,'--','LineWidth',linewidth,'Color',[139,69,19]/256);
hold off;  
xlabel('Offshore distance (km)');
ylabel('Elevation (km)');
set(gca,'FontSize',fontsize);
axis([0 250 -4 0]);



cntr = cntr + 1;
ax = axes('position',axpos(cntr,:)); 
defaultcolororder = get(gca,'Colororder');
plot(rhoRef,zz/1000,'LineWidth',linewidth);
hold on;
plot(rhoProf_wind,zz/1000,'.');
for k=1:length(E0_wind)
  plot([min(rhoRef) max(rhoRef)],[E0_wind(k) E0_wind(k)]/1000,'--','Color',defaultcolororder(2,:));
end
hold off;
xlabel('Potential density (kg/m^3)');
set(gca,'FontSize',fontsize);
text(1024.1,-3.8,'Wind-forced experiments','FontSize',fontsize);

cntr = cntr + 1;
ax = axes('position',axpos(cntr,:));  
plot(rhoRef,zz/1000,'LineWidth',linewidth);
hold on;
plot(rhoProf_rand,zz/1000,'.');
for k=1:length(E0_rand)
  plot([min(rhoRef) max(rhoRef)],[E0_rand(k) E0_rand(k)]/1000,'--','Color',defaultcolororder(2,:));
end
hold off;
xlabel('Potential density (kg/m^3)');
set(gca,'FontSize',fontsize);
text(1024.1,-3.8,'Randomly forced experiments','FontSize',fontsize);


%%% Add axis labels
for cntr = 1:size(axpos,1)
  annotation('textbox',[axpos(cntr,1)-0.05 axpos(cntr,2)-0.05 lab_size],'String',axlabels{cntr},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
end


%%% Add annotations
figure1 = gcf;
% Create doublearrow
annotation(figure1,'doublearrow',[0.132380952380953 0.174444444444444],...
  [0.29580300335902 0.296666666666667]);

% Create doublearrow
annotation(figure1,'doublearrow',[0.117904337681621 0.165924276169265],...
  [0.20146338071751 0.201997780244173]);

% Create textbox
annotation(figure1,'textbox',...
  [0.167936507936509 0.280662469627847 0.0721428571428571 0.0385714285714284],...
  'String','W_s_h_e_l_f',...
  'LineStyle','none',...
  'FontSize',14,...
  'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
  [0.120052320854103 0.156356143323741 0.072142857142857 0.0385714285714284],...
  'String','W_c_a_n_y_o_n_s',...
  'LineStyle','none',...
  'FontSize',14,...
  'FitBoxToText','off');

% Create textbox
annotation(figure1,'textbox',...
  [0.209139180542299 0.192982114466915 0.072142857142857 0.0385714285714284],...
  'String','S_b_o_t',...
  'LineStyle','none',...
  'FontSize',14,...
  'FitBoxToText','off');

% Create line
annotation(figure1,'line',[0.161469933184855 0.207126948775056],...
  [0.257491675915649 0.220865704772475]);

% Create textbox
annotation(figure1,'textbox',...
  [0.117402198889949 0.737532309450308 0.0592857142857143 0.0385714285714287],...
  'String',{'\tau_m_a_x'},...
  'LineStyle','none',...
  'FontSize',14,...
  'FitBoxToText','off');

% Create doublearrow
annotation(figure1,'doublearrow',[0.0841542051118892 0.183360554318238],...
  [0.776068291665944 0.776020870207735]);

% Create doublearrow
annotation(figure1,'doublearrow',[0.187717679499417 0.286924028705765],...
  [0.776068291665944 0.776020870207735]);

% Create textbox
annotation(figure1,'textbox',...
  [0.210945876197546 0.774168146174936 0.0592857142857142 0.0385714285714286],...
  'String','P_x^0',...
  'LineStyle','none',...
  'FontSize',14,...
  'FitBoxToText','off');

% % Create doublearrow
% annotation(figure1,'doublearrow',[0.0841542051118892 0.183360554318238],...
%   [0.794957180554833 0.794909759096624]);
% 
% % Create doublearrow
% annotation(figure1,'doublearrow',[0.187717679499417 0.286924028705765],...
%   [0.794957180554833 0.794909759096624]);
% 
% % Create textbox
% annotation(figure1,'textbox',...
%   [0.116291087778838 0.756421198339197 0.0592857142857143 0.0385714285714286],...
%   'String',{'\tau_m_a_x'},...
%   'LineStyle','none',...
%   'FontSize',14,...
%   'FitBoxToText','off');
% 
% % Create textbox
% annotation(figure1,'textbox',...
%   [0.210945876197546 0.784168146174936 0.0592857142857142 0.0385714285714286],...
%   'String','F_b_a_r_o',...
%   'LineStyle','none',...
%   'FontSize',14,...
%   'FitBoxToText','off');
