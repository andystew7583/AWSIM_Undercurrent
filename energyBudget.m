%%%
%%% energyBudget.m
%%%
%%% Computes terms from the thickness-weighted energy budget. The energy
%%% conversion terms are split based on the gridpoints on which they are
%%% most naturally defined.
%%%

function energyBudget (local_home_dir, run_name)

  %%% Load simulation parameters and time-averaged output
  loadParams;
  load([run_name,'_tavg.mat']);
  
  %%% Calculate Gsum, the sum of reduced gravities down to the layer in
  %%% question. This definition gives gsum_k*rho0 = g*rho_k, so this is
  %%% effectively a measure of the density in each layer.
  Gsum = repmat(reshape(cumsum(gg),[1 1 Nlay]),[Nx Ny 1]);

  %%% Thickness-weighted average velocities
  uu_twa = hu_tavg ./ hh_w_tavg;
  vv_twa = hv_tavg ./ hh_s_tavg;
  usq_twa = husq_tavg ./ hh_w_tavg;
  vsq_twa = hvsq_tavg ./ hh_s_tavg;
  uu_bol = uu_twa - uu_tavg;
  vv_bol = vv_twa - vv_tavg;
  
  %%% Total/mean/eddy kinetic/potential energies
  TKE_w = 0.5.*hh_w_tavg.*usq_twa;
  TKE_s = 0.5*hh_s_tavg.*vsq_twa;
  TKE = 0.5 * (TKE_w(1:Nx,:,:) + TKE_w([2:Nx 1],:,:)) + 0.5 * (TKE_s(:,1:Ny,:) + TKE_s(:,[2:Ny 1],:));
  MKE_w = 0.5*hh_w_tavg.*uu_twa.^2; 
  MKE_s = 0.5*hh_s_tavg.*vv_twa.^2;
  MKE = 0.5 * (MKE_w(1:Nx,:,:) + MKE_w([2:Nx 1],:,:)) + 0.5 * (MKE_s(:,1:Ny,:) + MKE_s(:,[2:Ny 1],:));
  EKE_w = TKE_w - MKE_w;
  EKE_s = TKE_s - MKE_s;
  EKE = 0.5 * (EKE_w(1:Nx,:,:) + EKE_w([2:Nx 1],:,:)) + 0.5 * (EKE_s(:,1:Ny,:) + EKE_s(:,[2:Ny 1],:));  
  TPE = Gsum.*hz_tavg;
  MPE = Gsum.*hh_tavg.*zz_tavg;
  EPE = TPE - MPE;
  
  %%% MPE to MKE conversion terms
  dMdx_tavg = (MM_tavg(1:Nx,:,:)-MM_tavg([Nx 1:Nx-1],:,:))/dx;
  dMdy_tavg = (MM_tavg(:,1:Ny,:)-MM_tavg(:,[Ny 1:Ny-1],:))/dy;   
  dzdx_tavg = (zz_tavg(1:Nx,:,:)-zz_tavg([Nx 1:Nx-1],:,:))/dx;
  dzdy_tavg = (zz_tavg(:,1:Ny,:)-zz_tavg(:,[Ny 1:Ny-1],:))/dy;     
  MPE_MKE_I = - Gsum.*hh_tavg.*dzdt_tavg;
  MPE_MKE_II = - Gsum.*hh_w_tavg.*uu_tavg.*dzdx_tavg;   
  MPE_MKE_III = - Gsum.*hh_s_tavg.*vv_tavg.*dzdy_tavg;
  MPE_MKE_IV = - Gsum.*hh_w_tavg.*uu_bol.*dzdx_tavg;   
  MPE_MKE_V = - Gsum.*hh_s_tavg.*vv_bol.*dzdy_tavg; 
  MPE_MKE_tot = MPE_MKE_I+MPE_MKE_II+MPE_MKE_III+MPE_MKE_IV+MPE_MKE_V;
  MPE_MKE_Aiki_I = - hh_w_tavg.*uu_tavg.*dMdx_tavg;   
  MPE_MKE_Aiki_II = - hh_s_tavg.*vv_tavg.*dMdy_tavg;
  MPE_MKE_Aiki_III = - hh_w_tavg.*uu_bol.*dMdx_tavg;   
  MPE_MKE_Aiki_IV = - hh_s_tavg.*vv_bol.*dMdy_tavg; 
  MPE_MKE_Aiki_tot = MPE_MKE_Aiki_I+MPE_MKE_Aiki_II+MPE_MKE_Aiki_III+MPE_MKE_Aiki_IV;
  
  %%% EPE to EKE conversion terms
  EPE_EKE_I = - Gsum.*hdzdt_tavg - MPE_MKE_I;
  EPE_EKE_II = - Gsum.*hudzdx_tavg - MPE_MKE_II - MPE_MKE_IV;
  EPE_EKE_III = - Gsum.*hvdzdy_tavg - MPE_MKE_III - MPE_MKE_V;
  EPE_EKE_tot = EPE_EKE_I + EPE_EKE_II + EPE_EKE_III;
  EPE_EKE_Aiki_I = - hudMdx_tavg - MPE_MKE_Aiki_I - MPE_MKE_Aiki_III;
  EPE_EKE_Aiki_II = - hvdMdy_tavg - MPE_MKE_Aiki_II - MPE_MKE_Aiki_IV;
  EPE_EKE_Aiki_III = - uu_tavg .* (hdMdx_tavg - hh_w_tavg.*dMdx_tavg);
  EPE_EKE_Aiki_IV = - vv_tavg .* (hdMdy_tavg - hh_s_tavg.*dMdy_tavg);
  EPE_EKE_Aiki_I = EPE_EKE_Aiki_I - EPE_EKE_Aiki_III;
  EPE_EKE_Aiki_II = EPE_EKE_Aiki_II - EPE_EKE_Aiki_IV;
  EPE_EKE_Aiki_tot = EPE_EKE_Aiki_I + EPE_EKE_Aiki_II + EPE_EKE_Aiki_III + EPE_EKE_Aiki_IV;
  
  %%% MKE to EKE conversion terms
  hphi_mean = hh_tavg.*phi_tavg;
  hphi_eddy = hphi_tavg - hphi_mean;
  pdedx_mean = 0.5*(phi_int_tavg(1:Nx,:,:)+phi_int_tavg([Nx 1:Nx-1],:,:)).*(eta_tavg(1:Nx,:,:)-eta_tavg([Nx 1:Nx-1],:,:))/dx;
  pdedx_eddy = pdedx_tavg - pdedx_mean;
  pdedy_mean = 0.5*(phi_int_tavg(:,1:Ny,:)+phi_int_tavg(:,[Ny 1:Ny-1],:)).*(eta_tavg(:,1:Ny,:)-eta_tavg(:,[Ny 1:Ny-1],:))/dy;
  pdedy_eddy = pdedy_tavg - pdedy_mean;
  MKE_EKE_I = uu_twa .* ((hphi_eddy(1:Nx,:,:)-hphi_eddy([Nx 1:Nx-1],:,:))/dx + diff(pdedx_eddy,1,3)); 
  MKE_EKE_II = vv_twa .* ((hphi_eddy(:,1:Ny,:)-hphi_eddy(:,[Ny 1:Ny-1],:))/dy + diff(pdedy_eddy,1,3)); 
  husq_eddy = husq_tavg - hh_w_tavg.*uu_twa.^2;
  hvsq_eddy = hvsq_tavg - hh_s_tavg.*vv_twa.^2;
  huv_eddy = huv_tavg - 0.5.*(uu_twa(1:Nx,:,:)+uu_twa([Nx 1:Nx-1],:,:)) ...
                .* 0.5.*(vv_twa(:,1:Ny,:)+vv_twa(:,[Ny 1:Ny-1],:)) ...
                .* 0.25.*(hh_tavg(1:Nx,1:Ny,:)+hh_tavg(1:Nx,[Ny 1:Ny-1],:)+hh_tavg([Nx 1:Nx-1],1:Ny,:)+hh_tavg([Nx 1:Nx-1],[Ny 1:Ny-1],:));                  
  MKE_EKE_III = - 0.5*(husq_eddy([2:Nx 1],:,:)+husq_eddy(1:Nx,:,:)).*(uu_twa([2:Nx 1],:,:)-uu_twa(1:Nx,:,:))/dx ...
                - 0.5*(hvsq_eddy(:,[2:Ny 1],:)+hvsq_eddy(:,1:Ny,:)).*(vv_twa(:,[2:Ny 1],:)-vv_twa(:,1:Ny,:))/dy;
  MKE_EKE_IV = - huv_eddy .* ( (uu_twa(:,1:Ny,:)-uu_twa(:,[Ny 1:Ny-1],:))/dy ...
                             + (vv_twa(1:Nx,:,:)-vv_twa([Nx 1:Nx-1],:,:))/dx );  
  MKE_EKE_tot = MKE_EKE_I + MKE_EKE_II + MKE_EKE_III + MKE_EKE_IV;
  MKE_EKE_Aiki_I = uu_twa .* (hdMdx_tavg - hh_w_tavg.*dMdx_tavg);
  MKE_EKE_Aiki_II = vv_twa .* (hdMdy_tavg - hh_s_tavg.*dMdy_tavg);
  MKE_EKE_Aiki_tot = MKE_EKE_Aiki_I + MKE_EKE_Aiki_II + MKE_EKE_III + MKE_EKE_IV;
  
  %%% Store for later
  save([run_name,'_energy.mat'], ...
    'uu_twa','vv_twa','usq_twa','vsq_twa','uu_bol','vv_bol', ...
    'TKE','MKE','EKE','TPE','MPE','EPE', ...
    'MPE_MKE_I','MPE_MKE_II','MPE_MKE_III','MPE_MKE_IV','MPE_MKE_V','MPE_MKE_tot', ...
    'MPE_MKE_Aiki_I','MPE_MKE_Aiki_II','MPE_MKE_Aiki_III','MPE_MKE_Aiki_IV','MPE_MKE_Aiki_tot', ...
    'EPE_EKE_I','EPE_EKE_II','EPE_EKE_III','EPE_EKE_tot', ...
    'EPE_EKE_Aiki_I','EPE_EKE_Aiki_II','EPE_EKE_Aiki_III','EPE_EKE_Aiki_IV','EPE_EKE_Aiki_tot', ...
    'MKE_EKE_I','MKE_EKE_II','MKE_EKE_III','MKE_EKE_IV','MKE_EKE_tot', ...
    'MKE_EKE_Aiki_I','MKE_EKE_Aiki_II','MKE_EKE_Aiki_tot');
  
end