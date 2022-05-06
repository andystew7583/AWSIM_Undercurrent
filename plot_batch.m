%%%
%%% plot_batch.m
%%%
%%% Plots diagnostics from our batch of simulations.
%%%

%%% Load constant parameters
constants;

%%% Location of runs on file system
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';

%%% Load parameters for all simulations in the batch
[is_wind_batch,grid_size_batch,wind_stress_batch, ...
  rand_force_batch,num_canyons_batch,amp_canyons_batch, ...
  max_slope_batch,sb_width_batch,baro_force_batch, ...
  drag_coeff_batch,end_time_batch] = loadBatchParams ('batchData_hires.txt');

%%% Iterate through simulations 
Nsims = length(is_wind_batch);
uc_speed_batch = zeros(Nsims,1);
uc_EKE_batch = zeros(Nsims,1);
tot_EKE_batch = zeros(Nsims,1);
uc_advection_batch = zeros(Nsims,1);
uc_eddyforce_batch = zeros(Nsims,1);
uc_hypervisc_batch = zeros(Nsims,1);
tot_advection_batch = zeros(Nsims,1);
uc_windStress_batch = zeros(Nsims,1);
uc_baroForcing_batch = zeros(Nsims,1);
uc_quadBotDrag_batch = zeros(Nsims,1);
uc_TFS_batch = zeros(Nsims,1);
for n=1:Nsims
  
  %%% Derived parameters
  if (is_wind_batch(n))
    config = 'wind';
    uc_layidx = 3;
  else
    config = 'rand';
    uc_layidx = 2;
  end
  grid_size = grid_size_batch(n);
  wind_stress = wind_stress_batch(n);
  rand_force = rand_force_batch(n);
  num_canyons = num_canyons_batch(n);
  amp_canyons = amp_canyons_batch(n);
  max_slope = max_slope_batch(n);
  sb_width = sb_width_batch(n);
  baro_force = baro_force_batch(n);
  drag_coeff = drag_coeff_batch(n);
  end_time = end_time_batch(n);
  
  %%% Generate simulation name
  run_name = constructRunName (config,grid_size,wind_stress, ...
          rand_force,num_canyons,amp_canyons,max_slope,sb_width,baro_force,drag_coeff);
  disp(['Working on ',run_name,' ...']);
        
  %%% Calculate along-isobath mean flow 
  loadParams;
  load(fullfile('products',[run_name,'_meanFlow.mat']));
  load(fullfile('products',[run_name,'_momBalance.mat']));
  
  
    
  
  %%% Compute diagnostics
  uc_didx = find((dd>=150) & (dd<=1000));
%   uc_didx = find((dd>=120) & (dd<=400));
  uc_yidx = find((yy_h>50*m1km) & (yy_h<100*m1km));
  uc_speed_batch(n) = avg_trap(y_avg(uc_didx),circ_twa(uc_didx,uc_layidx));
%   uc_EKE_batch(n) = mean(EKE(uc_yidx,uc_layidx)) / mean(h_avg(uc_yidx,uc_layidx));
  uc_EKE_batch(n) = avg_trap(y_avg(uc_didx),sum(EKE(uc_didx,uc_layidx:end),2)) / avg_trap(y_avg(uc_didx),sum(h_avg(uc_didx,uc_layidx:end),2));
  tot_EKE_batch(n) = avg_trap(y_avg(uc_didx),sum(EKE(uc_didx,1:end),2)) / avg_trap(y_avg(uc_didx),sum(h_avg(uc_didx,1:end),2));
  uc_advection_batch(n) = avg_trap(y_avg(uc_didx),sum(circ_advection(uc_didx,uc_layidx:end),2)./cntrlen(uc_didx));
  tot_advection_batch(n) = avg_trap(y_avg(uc_didx),sum(circ_advection(uc_didx,1:end),2)./cntrlen(uc_didx));
  uc_quadBotDrag_batch(n) = avg_trap(y_avg(uc_didx),sum(circ_quadBotDrag(uc_didx,uc_layidx:end),2)./cntrlen(uc_didx));
  uc_baroForcing_batch(n) = avg_trap(y_avg(uc_didx),sum(circ_baroForcing(uc_didx,uc_layidx:end),2)./cntrlen(uc_didx));
  uc_windStress_batch(n) = avg_trap(y_avg(uc_didx),sum(circ_windStress(uc_didx,uc_layidx:end),2)./cntrlen(uc_didx));
  uc_hypervisc_batch(n) = avg_trap(y_avg(uc_didx),sum(circ_hypervisc(uc_didx,uc_layidx:end),2)./cntrlen(uc_didx));
  uc_eddyforce_batch(n) = avg_trap(y_avg(uc_didx),sum(circ_advection(uc_didx,uc_layidx:end),2)./cntrlen(uc_didx)+sum(circ_Montgomery(uc_didx,uc_layidx:end),2)./cntrlen(uc_didx));
  uc_TFS_batch(n) = mean(mean(sum(UMom_Montgomery(:,uc_yidx,:),3),1));
end


figure(51);
scatter(uc_EKE_batch,uc_speed_batch,20,'b')

figure(52);
scatter(uc_advection_batch,uc_speed_batch,20,'b')

figure(53);
idx = find(is_wind_batch==1 & baro_force_batch==0);
scatter(uc_TFS_batch(idx).^.5,uc_speed_batch(idx),20,'b')
hold on;
idx = find(is_wind_batch==0 & baro_force_batch==0);
scatter(uc_TFS_batch(idx).^.5,uc_speed_batch(idx),20,'r')
hold off

figure(54);
idx = find(is_wind_batch==1 & baro_force_batch==0);
scatter(0.6*uc_EKE_batch(idx).^.5,uc_speed_batch(idx),20,'b')
% scatter(0.7*uc_EKE_batch(idx).^.5.*(5./sb_width_batch(idx)),uc_speed_batch(idx),20,'b')
hold on;
idx = find(is_wind_batch==0 & baro_force_batch==0);
scatter(0.6*uc_EKE_batch(idx).^.5,uc_speed_batch(idx),20,'r')
% scatter(0.7*uc_EKE_batch(idx).^.5.*(5./sb_width_batch(idx)),uc_speed_batch(idx),20,'r')
idx = find(baro_force_batch~=0);
scatter(0.5*uc_EKE_batch(idx).^.5,uc_speed_batch(idx),20,'g')
plot([0 0.1],[0 0.1],'k--');
hold off


figure(56);
idx = find(is_wind_batch==1 & baro_force_batch==0);
scatter(tot_EKE_batch(idx)./sb_width_batch(idx),tot_advection_batch(idx),20,'b');
hold on;
idx = find(is_wind_batch==0 & baro_force_batch==0);
scatter(tot_EKE_batch(idx)./sb_width_batch(idx),tot_advection_batch(idx),20,'r');
hold off

figure(57);
idx = find(is_wind_batch==1 & baro_force_batch==0);
scatter(uc_EKE_batch(idx)./sb_width_batch(idx),uc_advection_batch(idx),20,'b');
hold on;
idx = find(is_wind_batch==0 & baro_force_batch==0);
scatter(uc_EKE_batch(idx)./sb_width_batch(idx),uc_advection_batch(idx),20,'r');
idx = find(baro_force_batch ~= 0);
scatter(uc_EKE_batch(idx)./sb_width_batch(idx),uc_advection_batch(idx),20,'g');
hold off

figure(58);
idx = find(is_wind_batch==1 & baro_force_batch==0);
scatter(tot_EKE_batch(idx),uc_speed_batch(idx),20,'b')
hold on;
idx = find(is_wind_batch==0 & baro_force_batch==0);
scatter(tot_EKE_batch(idx),uc_speed_batch(idx),20,'r')
hold off


%%% Test eddy forcing scaling
figure(58);
eddyforce_theory = uc_EKE_batch./sb_width_batch;
eddyforce_theory_coeff = eddyforce_theory \ (uc_eddyforce_batch);
idx = find(is_wind_batch==1 & baro_force_batch==0);
scatter(eddyforce_theory_coeff.*eddyforce_theory(idx),uc_eddyforce_batch(idx),20,'b');
hold on;
idx = find(is_wind_batch==0 & baro_force_batch==0);
scatter(eddyforce_theory_coeff.*eddyforce_theory(idx),uc_eddyforce_batch(idx),20,'r');
hold off

%%% Test drag scaling: drag velocity is basically just sqrt(EKE), drag acts
%%% on undercurrent mean flow
drag_theory = drag_coeff_batch.*uc_speed_batch.*uc_EKE_batch.^.5;
figure(60)
idx = 1:Nsims;
drag_theory_coeff = drag_theory(idx) \ (-uc_quadBotDrag_batch(idx));
scatter(drag_theory_coeff.*drag_theory(idx),-uc_quadBotDrag_batch(idx));

%%% Test mom balance: eddy forcing vs drag
figure(61)
idx = 1:Nsims;
scatter(-uc_quadBotDrag_batch(idx)-uc_hypervisc_batch(idx),uc_advection_batch(idx)+uc_baroForcing_batch(idx))

%%% Test mom balance: eddy forcing vs drag
figure(62)
idx = 1:Nsims;
scatter(-uc_quadBotDrag_batch(idx)-uc_hypervisc_batch(idx),uc_eddyforce_batch(idx)+uc_baroForcing_batch(idx)+uc_windStress_batch(idx),'b')
hold on;
idx = find(is_wind_batch==0 & baro_force_batch==0);
scatter(-uc_quadBotDrag_batch(idx)-uc_hypervisc_batch(idx),uc_eddyforce_batch(idx)+uc_baroForcing_batch(idx)+uc_windStress_batch(idx),'r')
idx = find(baro_force_batch ~= 0);
scatter(-uc_quadBotDrag_batch(idx)-uc_hypervisc_batch(idx),uc_eddyforce_batch(idx)+uc_baroForcing_batch(idx)+uc_windStress_batch(idx),'g')
hold off;

%%% Test mom balance: eddy forcing vs drag
figure(63)
idx = 1:Nsims;
scatter(-uc_quadBotDrag_batch(idx),uc_eddyforce_batch(idx)+uc_baroForcing_batch(idx),'b')
hold on;
idx = find(is_wind_batch==0 & baro_force_batch==0);
scatter(-uc_quadBotDrag_batch(idx),uc_eddyforce_batch(idx)+uc_baroForcing_batch(idx),'r')
idx = find(baro_force_batch ~= 0);
scatter(-uc_quadBotDrag_batch(idx),uc_eddyforce_batch(idx)+uc_baroForcing_batch(idx),'g')
hold off;

%%% Test theoretical momentum balance: eddy forcing vs bottom friction
figure(64);
idx = find(is_wind_batch==1 & baro_force_batch==0);
scatter(eddyforce_theory_coeff.*eddyforce_theory(idx),drag_theory_coeff.*drag_theory(idx),20,'b')
hold on;
idx = find(is_wind_batch==0 & baro_force_batch==0);
scatter(eddyforce_theory_coeff.*eddyforce_theory(idx),drag_theory_coeff.*drag_theory(idx),20,'r')
idx = find(baro_force_batch~=0);
scatter(eddyforce_theory_coeff.*eddyforce_theory(idx),drag_theory_coeff.*drag_theory(idx),20,'g')
hold off

%%% Test theoretical undercurrent speed, derived from theoretical momentum
%%% balance
uc_speed_theory = uc_EKE_batch.^.5./sb_width_batch./drag_coeff_batch.*eddyforce_theory_coeff./drag_theory_coeff;
figure(65);
idx = find(is_wind_batch==1 & baro_force_batch==0);
scatter(uc_speed_theory(idx),uc_speed_batch(idx),20,'b')
hold on;
idx = find(is_wind_batch==0 & baro_force_batch==0);
scatter(uc_speed_theory(idx),uc_speed_batch(idx),20,'r')
idx = find(baro_force_batch~=0);
scatter(uc_speed_theory(idx),uc_speed_batch(idx),20,'g')
plot([0 0.15],[0 0.15],'k--');
hold off


figure(66);
scatter(log10(eddyforce_theory_coeff.*eddyforce_theory),log10(drag_theory_coeff.*drag_theory),20,'g')


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

