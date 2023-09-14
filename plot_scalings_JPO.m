% %%%
% %%% plot_batch.m
% %%%
% %%% Plots scalings for eddy forcing of the undercurrent from our batch of simulations, for our JPO paper.
% %%%
% 
% %%% Load constant parameters
% constants;
% rho0 = 1000;
% 
% %%% Location of runs on file system
% local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';
% 
% %%% Load parameters for all simulations in the batch
% [is_wind_batch,grid_size_batch,wind_stress_batch, ...
%   rand_force_batch,num_canyons_batch,amp_canyons_batch, ...
%   max_slope_batch,sb_width_batch,baro_force_batch, ...
%   drag_coeff_batch,end_time_batch] = loadBatchParams ('batchData_hires.txt');
% 
% %%% Iterate through simulations 
% Nsims = length(is_wind_batch);
% uc_speed_batch = zeros(Nsims,1);
% uc_EKE_batch = zeros(Nsims,1);
% tot_EKE_batch = zeros(Nsims,1);
% uc_advection_batch = zeros(Nsims,1);
% uc_eddyforce_batch = zeros(Nsims,1);
% uc_hypervisc_batch = zeros(Nsims,1);
% tot_advection_batch = zeros(Nsims,1);
% uc_windStress_batch = zeros(Nsims,1);
% uc_baroForcing_batch = zeros(Nsims,1);
% uc_Montgomery_batch = zeros(Nsims,1);
% uc_quadBotDrag_batch = zeros(Nsims,1);
% uc_TFS_batch = zeros(Nsims,1);
% for n=1:Nsims
%   
%   %%% Derived parameters
%   if (is_wind_batch(n))
%     config = 'wind';
%   else
%     config = 'rand';
%   end  
%   grid_size = grid_size_batch(n);
%   wind_stress = wind_stress_batch(n);
%   rand_force = rand_force_batch(n);
%   num_canyons = num_canyons_batch(n);
%   amp_canyons = amp_canyons_batch(n);
%   max_slope = max_slope_batch(n);
%   sb_width = sb_width_batch(n);
%   baro_force = baro_force_batch(n);
%   drag_coeff = drag_coeff_batch(n);
%   end_time = end_time_batch(n);
%   
%   %%% Generate simulation name
%   run_name = constructRunName (config,grid_size,wind_stress, ...
%           rand_force,num_canyons,amp_canyons,max_slope,sb_width,baro_force,drag_coeff);
%   disp(['Working on ',run_name,' ...']);
%         
%   %%% Calculate along-isobath mean flow 
%   loadParams;
%   load(fullfile('products',[run_name,'_meanFlow.mat']));
%   load(fullfile('products',[run_name,'_momBalance.mat']));
%   
%   
%     
%   
%   %%% Range of depths (isobaths) over which to average  
% %   uc_didx = find((dd>=150) & (dd<=1000));
%   uc_didx = find((dd>=120) & (dd<=600));
%   
%   %%% Determine depth of the undercurrent based on max northward flow speed
%   ucmax = max(circ_twa(uc_didx,:),[],1);
%   uc_layidx = find(ucmax == max(ucmax));  
%   
%   %%% Range of latitudes over which to compute linear momentum budget
%   %%% diagnostics
%   uc_yidx = find((yy_h>(75-amp_canyons)*m1km) & (yy_h<(75+amp_canyons)*m1km));
% 
%   %%% Compute diagnostics  
%   uc_speed_batch(n) = avg_trap(y_avg(uc_didx),circ_twa(uc_didx,uc_layidx));
% %   uc_EKE_batch(n) = mean(EKE(uc_yidx,uc_layidx)) / mean(h_avg(uc_yidx,uc_layidx));
%   uc_EKE_batch(n) = avg_trap(y_avg(uc_didx),sum(EKE(uc_didx,uc_layidx:end).*h_avg(uc_didx,uc_layidx:end),2)) / avg_trap(y_avg(uc_didx),sum(h_avg(uc_didx,uc_layidx:end),2));
%   tot_EKE_batch(n) = avg_trap(y_avg(uc_didx),sum(EKE(uc_didx,1:end),2).*h_avg(uc_didx,uc_layidx:end)) / avg_trap(y_avg(uc_didx),sum(h_avg(uc_didx,1:end),2));
%   uc_advection_batch(n) = avg_trap(y_avg(uc_didx),sum(circ_advection(uc_didx,uc_layidx:end),2)./cntrlen(uc_didx));
%   tot_advection_batch(n) = avg_trap(y_avg(uc_didx),sum(circ_advection(uc_didx,1:end),2)./cntrlen(uc_didx));
%   uc_quadBotDrag_batch(n) = avg_trap(y_avg(uc_didx),sum(circ_quadBotDrag(uc_didx,uc_layidx:end),2)./cntrlen(uc_didx));
%   uc_baroForcing_batch(n) = avg_trap(y_avg(uc_didx),sum(circ_baroForcing(uc_didx,uc_layidx:end),2)./cntrlen(uc_didx));
%   uc_Montgomery_batch(n) = avg_trap(y_avg(uc_didx),sum(circ_Montgomery(uc_didx,uc_layidx:end),2)./cntrlen(uc_didx));
%   uc_windStress_batch(n) = avg_trap(y_avg(uc_didx),sum(circ_windStress(uc_didx,uc_layidx:end),2)./cntrlen(uc_didx));
%   uc_hypervisc_batch(n) = avg_trap(y_avg(uc_didx),sum(circ_hypervisc(uc_didx,uc_layidx:end),2)./cntrlen(uc_didx));
% %   uc_eddyforce_batch(n) = avg_trap(y_avg(uc_didx),sum(circ_advection(uc_didx,uc_layidx:end),2)./cntrlen(uc_didx)+sum(circ_Montgomery(uc_didx,uc_layidx:end),2)./cntrlen(uc_didx));
% %   uc_eddyforce_batch(n) = avg_trap(y_avg(uc_didx),sum(circ_advection(uc_didx,uc_layidx:end),2)./cntrlen(uc_didx));
%   uc_eddyforce_batch(n) = avg_trap(y_avg(uc_didx),sum(rho0*circ_eddyMomFlux(uc_didx,uc_layidx:end),2)./cntrlen(uc_didx));
%   uc_TFS_batch(n) = mean(mean(sum(UMom_Montgomery(:,uc_yidx,:),3),1));
% end
% 





eddyforce_theory = rho0.*uc_EKE_batch.*100./(sb_width_batch*1000);
% eddyforce_theory(wind_stress_batch>0) = eddyforce_theory(wind_stress_batch>0).*wind_stress_batch(wind_stress_batch>0)./0.05;
eddyforce_theory_coeff = eddyforce_theory \ (uc_eddyforce_batch);


% drag_theory = rho0.*drag_coeff_batch*1e-3.*uc_speed_batch.*sqrt(2*uc_EKE_batch+uc_speed_batch.^2);
drag_theory = rho0.*drag_coeff_batch*1e-3.*uc_speed_batch.*sqrt(2*uc_EKE_batch);
% drag_theory_coeff = drag_theory(idx) \ (-uc_quadBotDrag_batch(idx));
drag_theory_coeff = 1;

% uc_speed_theory = (uc_EKE_batch.^.5)./sb_width_batch./drag_coeff_batch.*eddyforce_theory_coeff./drag_theory_coeff;
% uc_speed_theory = (eddyforce_theory_coeff.*eddyforce_theory)./(drag_theory_coeff.*1e-3.*rho0.*drag_coeff_batch.*sqrt(2*uc_EKE_batch));
uc_speed_theory = (eddyforce_theory_coeff.*eddyforce_theory+uc_baroForcing_batch/1)./(drag_theory_coeff.*1e-3.*rho0.*drag_coeff_batch.*sqrt(2*uc_EKE_batch));


%%% Plotting options
fontsize = 14;
axpos = zeros(4,4);
axpos(1,:) = [0.075 0.57 .40 .41];
axpos(2,:) = [0.575 0.57 .40 .41];
axpos(3,:) = [0.075 0.07 .40 .41];
axpos(4,:) = [0.575 0.07 .40 .41];
axlabels = {'(a)','(b)','(c)','(d)'};
lab_size = [0.05 0.03];
rho0 = 1000;
markersize = 40;

defaultcolororder = get(gca,'ColorOrder');
% colororder = zeros(8,3);
% colororder(1:idx_ref_Cd-1,:) = defaultcolororder(1:idx_ref_Cd-1,:);
% colororder(idx_ref_Cd,:) = [0 0 0];
% colororder(idx_ref_Cd+1:N_Cd,:) = defaultcolororder(idx_ref_Cd:N_Cd-1,:);
% 
% markersize = 10;
% markershapes = {'>','o','*','<','v','d','^','s','x','+'};


figure(109);
clf;
set(gcf,'Position',[382   306   920   720]);

defaultcolororder = get(gca,'ColorOrder');

idx_wind_stress = find(is_wind_batch==1 & baro_force_batch==0 & num_canyons_batch==4 & amp_canyons_batch==25 & drag_coeff_batch==2 & max_slope_batch==0.15 & sb_width_batch==5);
idx_wind_canyons = find(is_wind_batch==1 & baro_force_batch==0 & num_canyons_batch~=4);
idx_wind_amp = find(is_wind_batch==1 & baro_force_batch==0 & amp_canyons_batch~=25);
idx_wind_drag = find(is_wind_batch==1 & baro_force_batch==0 & drag_coeff_batch~=2);
idx_wind_slope = find(is_wind_batch==1 & baro_force_batch==0 & max_slope_batch~=0.15);
idx_wind_width = find(is_wind_batch==1 & baro_force_batch==0 & sb_width_batch~=5);
idx_wind_baro = find(baro_force_batch~=0);
idx_rand_force = find(is_wind_batch==0 & num_canyons_batch==4 & amp_canyons_batch==25 & drag_coeff_batch==2 & max_slope_batch==0.15 & sb_width_batch==5);
idx_rand_canyons = find(is_wind_batch==0 & num_canyons_batch~=4);
idx_rand_amp = find(is_wind_batch==0 & amp_canyons_batch~=25);
idx_rand_drag = find(is_wind_batch==0 & drag_coeff_batch~=2);
idx_rand_slope = find(is_wind_batch==0 & max_slope_batch~=0.15);
idx_rand_width = find(is_wind_batch==0 & sb_width_batch~=5);

idx_set = {idx_wind_stress,idx_wind_canyons,idx_wind_amp,idx_wind_drag,idx_wind_slope,idx_wind_width,...
            idx_rand_force,idx_rand_canyons,idx_rand_amp,idx_rand_drag,idx_rand_slope,idx_rand_width,idx_wind_baro,[4],[25]};
facecolors = [repmat(defaultcolororder(1,:),[6 1]);repmat(defaultcolororder(2,:),[6 1]);defaultcolororder(3,:);[0 0 80/255];[80/255 0 0]];
markerstyles = ['o','^','v','s','<','d','o','^','v','s','<','d','p','o','o'];
legstrs={'$\tau_{\mathrm{max}}$','$N_{\mathrm{canyons}}$','$W_{\mathrm{canyons}}$','$C_{\mathrm{drag}}$',...
  '$S_{\mathrm{bot}}$','$W_{\mathrm{shelf}}$','$F_{\mathrm{rand}}$','$N_{\mathrm{canyons}}$','$W_{\mathrm{canyons}}$',...
  '$C_{\mathrm{drag}}$','$S_{\mathrm{bot}}$','$W_{\mathrm{shelf}}$','$P_x^0$','Ref. $\tau_{\mathrm{max}}$','Ref. $F_{\mathrm{rand}}$'};


%%% Test eddy forcing scaling
subplot('Position',axpos(1,:));
for i=1:length(idx_set)
  idx = idx_set{i};
  scatter(eddyforce_theory_coeff.*eddyforce_theory(idx),uc_eddyforce_batch(idx),markersize,markerstyles(i),'MarkerFaceColor',facecolors(i,:),'MarkerEdgeColor','k');
  if (i == 1)   
    hold on;
  end
end
plot([0 0.014],[0 0.014],'k--');
hold off
axis([0 0.014 0 0.014]);
xlabel('Eddy force, theory (N/m^2)');
ylabel('Eddy force, diagnosed (N/m^2)');
set(gca,'FontSize',14);
handle = legend(legstrs,'Location','SouthEast','interpreter','latex');
set(handle,'Position',[0.388889059813126 0.58125 0.1035022445347 0.330277775393592]);
r = corr(eddyforce_theory,uc_eddyforce_batch);
text(1e-3,1e-2,['r^2 = ',num2str(r^2)],'FontSize',14);

%%% Test drag scaling: drag velocity is basically just sqrt(EKE), drag acts
%%% on undercurrent mean flow
subplot('Position',axpos(2,:));
for i=1:length(idx_set)
  idx = idx_set{i};
  scatter(drag_theory_coeff.*drag_theory(idx),-uc_quadBotDrag_batch(idx),markersize,markerstyles(i),'MarkerFaceColor',facecolors(i,:),'MarkerEdgeColor','k');
  if (i == 1)   
    hold on;
  end
end
plot([0 0.014],[0 0.014],'k--');
hold off
axis([0 0.014 0 0.014]);
xlabel('Drag force, theory (N/m^2)');
ylabel('Drag force, diagnosed (N/m^2)');
set(gca,'FontSize',14);
r = corr(drag_theory,-uc_quadBotDrag_batch);
text(1e-3,1e-2,['r^2 = ',num2str(r^2)],'FontSize',14);

%%% Test theoretical momentum balance: eddy forcing vs bottom friction
subplot('Position',axpos(3,:));
for i=1:length(idx_set)
  idx = idx_set{i};
%   scatter(eddyforce_theory_coeff.*eddyforce_theory(idx),drag_theory_coeff.*drag_theory(idx),markersize,markerstyles(i),'MarkerFaceColor',facecolors(i,:),'MarkerEdgeColor','k');
  scatter(eddyforce_theory_coeff.*eddyforce_theory(idx)+uc_baroForcing_batch(idx)/1,drag_theory_coeff.*drag_theory(idx),markersize,markerstyles(i),'MarkerFaceColor',facecolors(i,:),'MarkerEdgeColor','k');
  if (i == 1)   
    hold on;
  end
end
plot([0 0.014],[0 0.014],'k--');
hold off
axis([0 0.014 0 0.014]);
xlabel('Eddy force, theory + APF (N/m^2)');
ylabel('Drag force, theory (N/m^2)');
set(gca,'FontSize',14);
% r = corr(eddyforce_theory_coeff.*eddyforce_theory,drag_theory_coeff.*drag_theory);
r = corr(eddyforce_theory_coeff.*eddyforce_theory+uc_baroForcing_batch/1,drag_theory_coeff.*drag_theory);
text(1e-3,1e-2,['r^2 = ',num2str(r^2)],'FontSize',14);

%%% Test theoretical undercurrent speed, derived from theoretical momentum
%%% balance
subplot('Position',axpos(4,:));
for i=1:length(idx_set)
  idx = idx_set{i};
  scatter(uc_speed_theory(idx),uc_speed_batch(idx),markersize,markerstyles(i),'MarkerFaceColor',facecolors(i,:),'MarkerEdgeColor','k');
  if (i == 1)   
    hold on;
  end
end
plot([0 0.1],[0 0.1],'k--');
hold off
xlabel('Undercurrent speed, theory (m/s)');
ylabel('Undercurrent speed, diagnosed (m/s)');
set(gca,'FontSize',14);
r = corr(uc_speed_theory,uc_speed_batch);
text(1e-2,.8e-1,['r^2 = ',num2str(r^2)],'FontSize',14);

%%% Add axis labels
for cntr = 1:size(axpos,1)
  annotation('textbox',[axpos(cntr,1)-0.05 axpos(cntr,2)-0.03 lab_size],'String',axlabels{cntr},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
end









%%% Plotting options
fontsize = 14;
axpos = zeros(2,4);
axpos(1,:) = [0.075 0.14 .40 .81];
axpos(2,:) = [0.575 0.14 .40 .81];
axlabels = {'(a)','(b)'};
lab_size = [0.05 0.03];
rho0 = 1000;


figure(110);
clf;
set(gcf,'Position',[382   306   920   400]);

subplot('Position',axpos(1,:));
for i=1:length(idx_set)
  idx = idx_set{i};
  scatter(uc_TFS_batch(idx).^1,uc_speed_batch(idx),markersize,markerstyles(i),'MarkerFaceColor',facecolors(i,:),'MarkerEdgeColor','k');
  if (i == 1)   
    hold on;
  end
end
idx_rand = [idx_rand_force; idx_rand_canyons; idx_rand_amp; idx_rand_drag; idx_rand_slope; idx_rand_width];
% Crand = uc_TFS_batch(idx_rand) \ uc_speed_batch(idx_rand);
Crand = 1 / (uc_speed_batch(idx_rand) \ uc_TFS_batch(idx_rand));
idx = [idx_wind_stress; idx_wind_canyons; idx_wind_amp; idx_wind_drag; idx_wind_slope; idx_wind_width];
Csq = uc_TFS_batch(idx) \ uc_speed_batch(idx).^2;
plot([0:0.001:0.04],sqrt(Csq)*[0:0.001:0.04].^0.5,':','Color',defaultcolororder(1,:),'LineWidth',2);
plot([0:0.001:0.01],Crand*[0:0.001:0.01],':','Color',defaultcolororder(2,:),'LineWidth',2);
plot([0 0],[-0.01 0.08],'k--');
plot([-0.02 0.04],[0 0],'k--');
hold off
text(0.03,0.055,'~TFS^0^.^5','FontSize',fontsize,'Color',defaultcolororder(1,:));
text(0.01,0.07,'r=0.53','FontSize',fontsize,'Color',defaultcolororder(2,:));
xlabel('Topographic form stress (N/m^2)');
ylabel('Undercurrent speed (m/s)');
handle = legend(legstrs,'Location','SouthWest','interpreter','latex');
set(handle,'Position',[0.0771739130434787 0.1425 0.103502244534699 0.594499995708466]);
set(gca,'FontSize',14);
set(gca,'YLim',[0 0.08]);

subplot('Position',axpos(2,:));
for i=[1:6 13 14]
  idx = idx_set{i};
  scatter(wind_stress_batch(idx),uc_EKE_batch(idx),markersize,markerstyles(i),'MarkerFaceColor',facecolors(i,:),'MarkerEdgeColor','k');
  if (i == 1)   
    hold on;
  end
end
C = wind_stress_batch(idx) \ uc_EKE_batch(idx);
plot([0 0.1],1.1*C*[0 0.1],'k:');
text(0.085,0.002,'~\tau_m_a_x^1','FontSize',fontsize);
% plot([-0.02 0.04],[0 0],'k--');
hold off
xlabel('Wind stress max. (N/m^2)');
ylabel('Undercurrent EKE (m^2/s^2)');
set(gca,'FontSize',14);


%%% Add axis labels
for cntr = 1:size(axpos,1)
  annotation('textbox',[axpos(cntr,1)-0.07 axpos(cntr,2)-0.05 lab_size],'String',axlabels{cntr},'interpreter','latex','FontSize',fontsize+2,'LineStyle','None');
end








%%%
%%% Convenience function to average on a non-uniform grid via the
%%% trapezoidal rule.
%%%
function f_avg = avg_trap (y,f)
 
  f_avg = int_trap(y,f) / int_trap(y,ones(size(f)));
 
end

%%%
%%% Convenience function to integrate on a non-uniform grid via the
%%% trapezoidal rule.
%%%
function f_int = int_trap (y,f)
 
  f_int = 0.5.*sum((f(2:end)+f(1:end-1)).*(y(2:end)-y(1:end-1))); 
 
end

