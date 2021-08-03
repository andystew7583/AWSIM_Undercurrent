%%%
%%% plotMomBudget.m
%%%
%%% Plots diagnostics from time-averaged momentum budget for our
%%% undercurrent runs.
%%%

%%% Select experiment
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';
run_name = 'undercurrent_RF_F0.05_Nc4_Ht200';
% run_name = 'undercurrent_WF_tau0.05_Nc4_localwind';
% run_name = 'undercurrent_WF_tau0.05_Nc4_Ht200';

%%% Load experiment parameters
dirpath = fullfile(local_home_dir,run_name);
loadParams;
tmin = 5.05*t1year;
tmax = 10.05*t1year;

%%% Define isobath grid
dd = [100:10:4000]';
Nd = length(dd);

%%% Average u-momentum diagnostics  
UMom_PVadvection = do_avg(dirpath,OUTN_UMOM_Q,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_Montgomery = do_avg(dirpath,OUTN_UMOM_GRADM,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_KEgradient = do_avg(dirpath,OUTN_UMOM_GRADKE,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_dhdt = do_avg(dirpath,OUTN_UMOM_DHDT,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_hypervisc = do_avg(dirpath,OUTN_UMOM_A4,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_quadBotDrag = do_avg(dirpath,OUTN_UMOM_CDBOT,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_windStress = do_avg(dirpath,OUTN_UMOM_WIND,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_relaxation = do_avg(dirpath,OUTN_UMOM_RELAX,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_baroForcing = do_avg(dirpath,OUTN_UMOM_FBARO,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_randomForcing = do_avg(dirpath,OUTN_UMOM_RAND,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_diapycnal = do_avg(dirpath,OUTN_UMOM_WDIA,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
% UMom_buoyForce = do_avg(dirpath,OUTN_UMOM_BUOY,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
% UMom_linBotDrag = do_avg(dirpath,OUTN_UMOM_RDRAG,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
% UMom_linSurfDrag = do_avg(dirpath,OUTN_UMOM_RSURF,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
% UMom_quadSurfDrag = do_avg(dirpath,OUTN_UMOM_CDSURF,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
UMom_adv = UMom_PVadvection + UMom_KEgradient - UMom_dhdt; %%% TODO fix!!!

%%% Average v-momentum diagnostics  
VMom_PVadvection = do_avg(dirpath,OUTN_VMOM_Q,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_Montgomery = do_avg(dirpath,OUTN_VMOM_GRADM,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_KEgradient = do_avg(dirpath,OUTN_VMOM_GRADKE,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_dhdt = do_avg(dirpath,OUTN_VMOM_DHDT,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_hypervisc = do_avg(dirpath,OUTN_VMOM_A4,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_quadBotDrag = do_avg(dirpath,OUTN_VMOM_CDBOT,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_windStress = do_avg(dirpath,OUTN_VMOM_WIND,Nx,Ny,Nlay,n0_avg_hv,N_avg_hv,dt_avg_hv,tmin,tmax,startTime);
VMom_relaxation = do_avg(dirpath,OUTN_VMOM_RELAX,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
VMom_baroForcing = do_avg(dirpath,OUTN_VMOM_FBARO,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
VMom_randomForcing = do_avg(dirpath,OUTN_VMOM_RAND,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
VMom_diapycnal = do_avg(dirpath,OUTN_VMOM_WDIA,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
% VMom_buoyForce = do_avg(dirpath,OUTN_VMOM_BUOY,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
% VMom_linBotDrag = do_avg(dirpath,OUTN_VMOM_RDRAG,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
% VMom_linSurfDrag = do_avg(dirpath,OUTN_VMOM_RSURF,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
% VMom_quadSurfDrag = do_avg(dirpath,OUTN_VMOM_CDSURF,Nx,Ny,Nlay,n0_avg_hu,N_avg_hu,dt_avg_hu,tmin,tmax,startTime);
VMom_adv = VMom_PVadvection + VMom_KEgradient - VMom_dhdt;

%%% Calculate vorticity budget quantities
curl_windStress = calc_curl(UMom_windStress,VMom_windStress,dx,dy);
curl_quadBotDrag = calc_curl(UMom_quadBotDrag,VMom_quadBotDrag,dx,dy);
curl_Montgomery = calc_curl(UMom_Montgomery,VMom_Montgomery,dx,dy);
curl_hypervisc = calc_curl(UMom_hypervisc,VMom_hypervisc,dx,dy);
curl_adv = calc_curl(UMom_adv,VMom_adv,dx,dy);

%%% Interpolate to q-points
etab_q = 0*hhb;
etab_q(:,1) = 0.5*(hhb(:,1)+hhb([Nx 1:Nx-1],1));
etab_q(:,2:Ny) = 0.25*(hhb(:,1:Ny-1)+hhb([Nx 1:Nx-1],1:Ny-1)+hhb(:,2:Ny)+hhb([Nx 1:Nx-1],2:Ny));

%%% Compute circulation tendencies along isobaths
circ_windStress = calc_iso_circ(curl_windStress,dx,dy,dd,etab_q);
circ_quadBotDrag = calc_iso_circ(curl_quadBotDrag,dx,dy,dd,etab_q);
circ_Montgomery = calc_iso_circ(curl_Montgomery,dx,dy,dd,etab_q);
circ_hypervisc = calc_iso_circ(curl_hypervisc,dx,dy,dd,etab_q);
circ_adv = calc_iso_circ(curl_adv,dx,dy,dd,etab_q);

%%% Calculate contour lengths
cntrlen = calc_iso_int(ones(Nx,Ny),dx,dy,dd,etab_q);



figure(11);
plot(yy_h,mean(sum(UMom_windStress,3),1));
hold on;
plot(yy_h,mean(sum(UMom_Montgomery,3),1));
plot(yy_h,mean(sum(UMom_adv,3),1));
plot(yy_h,mean(sum(UMom_quadBotDrag,3),1));
plot(yy_h,mean(sum(UMom_hypervisc,3),1));
plot(yy_h,mean(sum(UMom_relaxation,3),1));
plot(yy_h,mean(sum(UMom_randomForcing,3),1));
plot(yy_h,mean(sum(UMom_diapycnal,3),1));
plot(yy_h,mean(sum(UMom_windStress+UMom_Montgomery+UMom_adv+UMom_quadBotDrag+UMom_hypervisc+UMom_relaxation+UMom_randomForcing+UMom_diapycnal,3),1));
hold off
legend('Wind','TFS','Adv','Drag','Visc','Sponge','Rand','Wdia','Total');

figure(12);
plot(yy_h,mean(sum(UMom_diapycnal,3),1));
hold on;
plot(yy_h,mean(sum(UMom_windStress+UMom_Montgomery+UMom_adv+UMom_quadBotDrag+UMom_hypervisc+UMom_relaxation+UMom_randomForcing+UMom_diapycnal,3),1));
hold off;

figure(13);
plot(dd,sum(circ_windStress,2)./cntrlen);
hold on;
plot(dd,sum(circ_Montgomery,2)./cntrlen);
plot(dd,sum(circ_adv,2)./cntrlen);
plot(dd,sum(circ_quadBotDrag,2)./cntrlen);
plot(dd,sum(circ_hypervisc,2)./cntrlen);
hold off;
legend('Wind','TFS','Adv','Drag','Visc');


layidx = 2;

figure(14);
plot(dd,sum(circ_windStress(:,1:layidx-1),2)./cntrlen);
hold on;
plot(dd,sum(circ_Montgomery(:,1:layidx-1),2)./cntrlen);
plot(dd,sum(circ_adv(:,1:layidx-1),2)./cntrlen);
plot(dd,sum(circ_quadBotDrag(:,1:layidx-1),2)./cntrlen);
plot(dd,sum(circ_hypervisc(:,1:layidx-1),2)./cntrlen);
hold off;
legend('Wind','TFS','Adv','Drag','Visc');

figure(15);
plot(dd,sum(circ_windStress(:,layidx:Nlay),2)./cntrlen);
hold on;
plot(dd,sum(circ_Montgomery(:,layidx:Nlay),2)./cntrlen);
plot(dd,sum(circ_adv(:,layidx:Nlay),2)./cntrlen);
plot(dd,sum(circ_quadBotDrag(:,layidx:Nlay),2)./cntrlen);
plot(dd,sum(circ_hypervisc(:,layidx:Nlay),2)./cntrlen);
hold off;
legend('Wind','TFS','Adv','Drag','Visc');


% plot(yy_h,mean(sum(UMom_baroForcing+UMom_buoyForce+UMom_linBotDrag+UMom_linSurfDrag+UMom_quadSurfDrag,3),1));
% plot(yy_h,mean(sum(UMom_diapycnal,3),1));
