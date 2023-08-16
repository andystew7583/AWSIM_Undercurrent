%%%
%%% anim.m
%%%
%%% Reads in data from the output of 'AWSIM' and makes a movie of the
%%% layer output.
%%%
function M = anim (local_home_dir,run_name,var,layer,tmin,tmax) 

  %%% Load parameters   
  loadParams;
  dirpath = fullfile(local_home_dir,run_name);
  
  %%% Make a movie of the data
  clf;
  axes('FontSize',16);
  M = moviein(1);      
  Mcnt = 1;
  
  %%% At each time iteration...
  for n=n0:1:n0+Nframes-1   
    
    %%% Current simulation time    
    t = startTime + (n-n0)*dt_s;
    tt(Mcnt) = t;
    
    if ((tmin >= 0 && t < tmin) || (tmax >= 0 && t > tmax))
      continue;
    end
    n
    
    %%% Set up figure
    figure(3); 
        
    switch (var)
      
      %%% Contour plot of relative vorticity
      case 'z'                
        
        %%% Load u
        data_file = fullfile(dirpath,[OUTN_U,num2str(layer-1),'_n=',num2str(n),'.dat']);
        uu = readOutputFile(data_file,Nx,Ny);
        
        %%% Load v
        data_file = fullfile(dirpath,[OUTN_V,num2str(layer-1),'_n=',num2str(n),'.dat']);
        vv = readOutputFile(data_file,Nx,Ny);            
        
        %%% Load h
        data_file = fullfile(dirpath,[OUTN_H,num2str(layer-1),'_n=',num2str(n),'.dat']);
        hh = readOutputFile(data_file,Nx,Ny); 
                     
        %%% Calculate the relative vorticity
        zeta = zeros(Nx+1,Ny+1);            
        zeta(2:Nx,2:Ny) = (uu(1:Nx-1,1:Ny-1)-uu(1:Nx-1,2:Ny)) / dy + (vv(2:Nx,2:Ny)-vv(1:Nx-1,2:Ny)) / dx;                      
        if (~useWallNS)
          zeta(:,Ny+1) = zeta(:,1);
        end
        zeta(Nx+1,:) = zeta(1,:);
        Ro = zeta./(2*Omega_z);
        max(max(abs(Ro)))        

        %%% Make the plot
        pcolor(XX_q/1000,YY_q/1000,Ro);
        shading interp;
        hold on
        contour(XX_h/1000,YY_h/1000,hhb,[-3600:600:-600],'EdgeColor','k');
        hold off;
        colorbar;
        colormap(cmocean('balance'));
        caxis([-0.4 0.4]);        
        title(strcat(['Ro=\zeta/|f0| at t=',num2str(t/t1year,'%.2f'),' years']));        
        xlabel('x (km)');
        ylabel('y (km)');                   
        set(gca,'FontSize',14);
        
      %%% Contour plot of potential vorticity
      case 'pv'
        
        %%% Load u
        data_file = fullfile(dirpath,[OUTN_U,num2str(layer-1),'_n=',num2str(n),'.dat']);
        uu = readOutputFile(data_file,Nx,Ny);
        
        %%% Load v
        data_file = fullfile(dirpath,[OUTN_V,num2str(layer-1),'_n=',num2str(n),'.dat']);
        vv = readOutputFile(data_file,Nx,Ny);    
        
        %%% Load h
        data_file = fullfile(dirpath,[OUTN_H,num2str(layer-1),'_n=',num2str(n),'.dat']);
        hh = readOutputFile(data_file,Nx,Ny);      
        
        %%% Calculate the relative vorticity
        zeta = zeros(Nx+1,Ny+1);            
        zeta(2:Nx,2:Ny) = (uu(1:Nx-1,1:Ny-1)-uu(1:Nx-1,2:Ny)) / dy + (vv(2:Nx,2:Ny)-vv(1:Nx-1,2:Ny)) / dx;                                         
        
        %%% Calculate PV        
        ff = 2*Omega_z;        
        hh_q = NaN*ones(Nx+1,Ny+1);
        hh_q(2:Nx,2:Ny) = 0.25 * (hh(2:Nx,1:Ny-1)+hh(2:Nx,2:Ny)+hh(1:Nx-1,1:Ny-1)+hh(1:Nx-1,2:Ny));        
        pv = (ff+zeta)./hh_q;                   
                
        %%% Make the plot
%         pcolor(XX_q/1000,YY_q/1000,pv);
        pcolor(XX_q/1000,YY_q/1000,log10(abs(pv)));
        shading interp;
        colorbar;
        colormap(pmkmp(100,'Swtth'));
        title(strcat(['PV at t=',num2str(t/t1year,'%.2f'),' years']));        
        xlabel('x (km)');
        ylabel('y (km)');         
        
      %%% Plot layer surface height
      case 'e'
        
        %%% Calculate layer surface height
        eta = hhb;
        for k=Nlay:-1:layer

          %%% Load kth layer thickness
          data_file = fullfile(dirpath,[OUTN_H,num2str(k-1),'_n=',num2str(n),'.dat']);
          hh = readOutputFile(data_file,Nx,Ny);      
          
          %%% Add layer thickness to eta
          eta = eta + hh;
          
        end
        
        %%% Make the plot
        pcolor(XX_h/1000,YY_h/1000,eta);
        shading interp;
        colorbar;
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']));
        colormap jet;        

      %%% Plot planetary pv
      case 'qp'
        
        %%% Load h
        data_file = fullfile(dirpath,[OUTN_H,num2str(layer-1),'_n=',num2str(n),'.dat']);
        hh = readOutputFile(data_file,Nx,Ny);    
        
        %%% Calculate planetary PV        
        ff = 2 * 2*Omega_z;
        hh_q = 0*ff;
        hh_q(1:Nx,2:Ny) = 0.25 * (hh(1:Nx,1:Ny-1)+hh(1:Ny,2:Ny)+hh([Nx 1:Nx-1],1:Ny-1)+hh([Nx 1:Nx-1],2:Ny));
        hh_q(1:Nx,Ny+1) = 0.5 * (hh(1:Ny,Ny)+hh([Nx 1:Nx-1],Ny));
        hh_q(1:Nx,1) = 0.5 * (hh(1:Nx,1)+hh([Nx 1:Nx-1],1));
        qp = ff./hh_q;
    
        %%% Make the plot
        pcolor(XX_q/1000,YY_q/1000,qp);
        shading interp;
        colorbar;
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']));
        colormap jet;       

      %%% Plot layer thickness
      case 'h'
        
        %%% Load h
        data_file = fullfile(dirpath,[OUTN_H,num2str(layer-1),'_n=',num2str(n),'.dat']);
        hh = readOutputFile(data_file,Nx,Ny);     
    
        %%% Make the plot
        pcolor(XX_h/1000,YY_h/1000,hh);
        shading interp;
        colorbar;
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']));
        colormap jet; 
%         caxis([0 1000])
        
      %%% Plot tracer
      case 'b'
        
        %%% Load h
        data_file = fullfile(dirpath,[OUTN_B,num2str(layer-1),'_n=',num2str(n),'.dat']);
        bb = readOutputFile(data_file,Nx,Ny);   
    
        %%% Make the plot
        pcolor(XX_h/1000,YY_h/1000,bb);
        shading interp;
        colorbar;
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']));
        colormap jet;
       
      %%% Plot surface pressure
      case 'p'
        
        %%% Calculate layer surface height
        eta = zeros(Nx,Ny,Nlay+1);
        hh = zeros(Nx,Ny,Nlay);
        eta(:,:,Nlay+1) = hhb;        
        for k=Nlay:-1:1

          %%% Load kth layer thickness
          data_file = fullfile(dirpath,[OUTN_H,num2str(k-1),'_n=',num2str(n),'.dat']);
          hh(:,:,k) = readOutputFile(data_file,Nx,Ny);    
          
          %%% Add layer thickness to eta
          eta(:,:,k) = eta(:,:,k+1) + hh(:,:,k);
          
        end
        
        %%% Calculate uppermost layer presure
        if (useRL)        
          %%% Load surface pressure
          data_file = fullfile(dirpath,[OUTN_PI,'_n=',num2str(n),'.dat']);
          pi = readOutputFile(data_file,Nx,Ny); 
        else
          pi = gg(1)*eta(:,:,1);
        end
        
        %%% Calculate hydrostatic pressure
        for k=2:layer
          pi = pi + gg(k)*eta(:,:,k);          
        end                
        
        gsum = sum(gg(1:layer));                
        pi = pi - gsum/3 * h0^4./hh(:,:,layer).^3;
    
        %%% Make the plot
        pcolor(XX_h/1000,YY_h/1000,pi);
        hold on
        [C,h] = contour(XX_h/1000,YY_h/1000,hhb,[-4000:500:-1000 -750],'EdgeColor','k');        
        caxis([min(min(pi)) max(max(pi))]);
        hold off;
        shading interp;
        colorbar;
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']));
        colormap jet;        
%         caxis([-200 200]);
       

      %%% Plot zonal velocity
      case 'u'
        
        %%% Load u
        data_file = fullfile(dirpath,[OUTN_U,num2str(layer-1),'_n=',num2str(n),'.dat']);
        uu = readOutputFile(data_file,Nx,Ny);       
    
        %%% Make the plot
        pcolor(XX_u/1000,YY_u/1000,uu);        
        shading interp;
        hold on
        [C,h] = contour(XX_h/1000,YY_h/1000,hhb,[-3600:600:-600],'EdgeColor','k');
        clabel(C,h);
        hold off;
        colorbar;
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']));
        colormap redblue;        
        caxis([-max(max(abs(uu))) max(max(abs(uu)))]);
%         caxis([-.1 .1]);
        
      %%% Plot meridional velocity
      case 'v'
        
        %%% Load v
        data_file = fullfile(dirpath,[OUTN_V,num2str(layer-1),'_n=',num2str(n),'.dat']);
        vv = readOutputFile(data_file,Nx,Ny);      
    
        %%% Make the plot
        pcolor(XX_v/1000,YY_v/1000,vv);
        shading interp;
        colorbar;
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']));
        colormap jet;        
        caxis([-max(max(abs(vv))) max(max(abs(vv)))]);
        
        
      case 'vcorr'
        
        %%% Load v
        uu = zeros(Nx,Ny,Nlay);
        vv = zeros(Nx,Ny,Nlay);
        for k=1:Nlay
          data_file = fullfile(dirpath,[OUTN_U,num2str(k-1),'_n=',num2str(n),'.dat']);
          uu(:,:,k) = readOutputFile(data_file,Nx,Ny);      
          data_file = fullfile(dirpath,[OUTN_V,num2str(k-1),'_n=',num2str(n),'.dat']);
          vv(:,:,k) = readOutputFile(data_file,Nx,Ny);      
        end
        uabs = sqrt(uu.^2+vv.^2);
        usurf = vv(:,:,1);
        usub = vv(:,:,3);
        scatter(usurf(:),usub(:))
        xlabel('Surface velocity');
        ylabel('Subsurface velocity');
        
        r = corr(usurf(:),usub(:));
        disp(r^2);

      %%% Plot zonally-averaged zonal velocity
      case 'ua'
        
        %%% Load u
        uu = zeros(Nx,Ny,Nlay);
        for k=1:Nlay
          data_file = fullfile(dirpath,[OUTN_U,num2str(k-1),'_n=',num2str(n),'.dat']);
          uu(:,:,k) = readOutputFile(data_file,Nx,Ny);    
        end        
        
        %%% Make the plot
        plot(yy_u/1000,mean(uu(:,:,1),1));        
        hold on
        for k=2:Nlay
          plot(yy_u/1000,mean(uu(:,:,k),1));        
        end
        plot(yy_u/1000,0*yy_u,'k--');
        hold off
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']));        
        axis([0 Ly/1000 -.3 .3]);
        
        
      %%% Plot zonally-averaged zonal velocity and relative vorticity
      case 'uaz'
        
        %%% Load u
        data_file = fullfile(dirpath,[OUTN_U,num2str(layer-1),'_n=',num2str(n),'.dat']);
        uu = readOutputFile(data_file,Nx,Ny);
        
        %%% Load v
        data_file = fullfile(dirpath,[OUTN_V,num2str(layer-1),'_n=',num2str(n),'.dat']);
        vv = readOutputFile(data_file,Nx,Ny);                           
                     
        %%% Calculate the relative vorticity
        zeta = zeros(Nx+1,Ny+1);            
        zeta(2:Nx,2:Ny) = (uu(1:Nx-1,1:Ny-1)-uu(1:Nx-1,2:Ny)) / dy + (vv(2:Nx,2:Ny)-vv(1:Nx-1,2:Ny)) / dx;                      
        if (~useWallNS)
          zeta(:,Ny+1) = zeta(:,1);
        end
        zeta(Nx+1,:) = zeta(1,:);
        Ro = zeta./(2*Omega_z);
        max(max(abs(Ro)))      
        
        H = mean(-hhb,1);
        dH_dy = 0*H;
        dH_dy(2:end-1) = (H(3:end)-H(1:end-2)) / (2*dy);        
        u_pred = dH_dy ./ (H.^3);
        u_pred = u_pred * -0.8/max(abs(u_pred));
        
        beta = dH_dy./H;
        idx_maxbeta = find(abs(beta)==max(abs(beta)));

        clf
        
        sum(u_pred.^2)
        sum(mean(uu,1).^2)
        
        %%% Make the plot
        ax = axes;
        set(ax,'Position',[0.05 0.08 0.2 0.85]);
        plot(mean(uu,1),yy_u/1000);        
        hold on        
        plot(u_pred,yy_u/1000,'k--');                
        plot(beta*-0.6/max(beta),yy_u/1000,'r--');        
        plot(0*yy_u,yy_u/1000,'k:');
%         plot([-1 1],[1 1]*yy_u(idx_maxbeta)/1000,'r--');
%         text(0.2,50,'$\beta_{\mathrm{max}}$','interpreter','latex','FontSize',14,'Color','r');
        hold off
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']),'interpreter','latex');        
        axis([-.6 .6 0 Ly/1000]);
        ylabel('y (km)','interpreter','latex');          
        xlabel('$\overline{u}^x$ (m/s)','interpreter','latex');
        set(gca,'FontSize',14);
        legend('Diagnosed','Theory','$\sim\beta_t$','interpreter','latex');
        
        %%% Make the plot
        ax = axes;
        set(ax,'Position',[0.32 0.08 0.65 0.85]);
        pcolor(XX_q/1000,YY_q/1000,Ro);
        shading interp;
        hold on
        [C,h] = contour(XX_h/1000,YY_h/1000,hhb,[-3800:600:-200],'EdgeColor','k');
        clabel(C,h);
        hold off;
        colorbar;
        colormap(cmocean('balance'));
        caxis([-0.4 0.4]);        
        title(strcat(['$Ro=\zeta/|f_0|$ at $t=$',num2str(t/t1day,'%.2f'),' days']),'interpreter','latex');        
        xlabel('$x$ (km)','interpreter','latex');
        ylabel('$y$ (km)','interpreter','latex');                   
        set(gca,'FontSize',14);
        
        
        %%% Plot zonally-averaged zonal velocity and relative vorticity
      case 'KEaz'
        
        %%% Load u
        data_file = fullfile(dirpath,[OUTN_U,num2str(layer-1),'_n=',num2str(n),'.dat']);
        uu = readOutputFile(data_file,Nx,Ny);
        
        %%% Load v
        data_file = fullfile(dirpath,[OUTN_V,num2str(layer-1),'_n=',num2str(n),'.dat']);
        vv = readOutputFile(data_file,Nx,Ny);            
        
        %%% Load h
        data_file = fullfile(dirpath,[OUTN_H,num2str(layer-1),'_n=',num2str(n),'.dat']);
        hh = readOutputFile(data_file,Nx,Ny);  
                     
        %%% Calculate the relative vorticity
        zeta = zeros(Nx+1,Ny+1);            
        zeta(2:Nx,2:Ny) = (uu(1:Nx-1,1:Ny-1)-uu(1:Nx-1,2:Ny)) / dy + (vv(2:Nx,2:Ny)-vv(1:Nx-1,2:Ny)) / dx;                      
        if (~useWallNS)
          zeta(:,Ny+1) = zeta(:,1);
        end
        zeta(Nx+1,:) = zeta(1,:);
        Ro = zeta./(2*Omega_z);
        max(max(abs(Ro)))      
        
        H = mean(-hhb,1);
        dH_dy = 0*H;
        dH_dy(2:end-1) = (H(3:end)-H(1:end-2)) / (2*dy);        
        KE_pred = dH_dy.^2 ./ (H.^5);
        KE_pred = KE_pred * 40/max(abs(KE_pred));

        clf
        
%         sum(u_pred.^2)
%         sum(mean(uu,1).^2)
        
        %%% Make the plot
        ax = axes;
        set(ax,'Position',[0.05 0.08 0.2 0.85]);
        plot(mean(0.5*hh.*(uu.^2+vv.^2),1),yy_u/1000);        
        hold on        
        plot(KE_pred,yy_u/1000,'k--');        
        plot(0*yy_u,yy_u/1000,'k:');
        hold off
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']),'interpreter','latex');        
%         axis([-.3 .3 0 Ly/1000]);
        ylabel('y (km)','interpreter','latex');          
        xlabel('$\overline{\frac{1}{2}h\mathbf{u}^2}^x$ (m/s)','interpreter','latex');
        set(gca,'FontSize',14);
        legend('Diagnosed','Theory');
        
        %%% Make the plot
        ax = axes;
        set(ax,'Position',[0.32 0.08 0.65 0.85]);
        pcolor(XX_q/1000,YY_q/1000,Ro);
        shading interp;
        hold on
        [C,h] = contour(XX_h/1000,YY_h/1000,hhb,[-3800:600:-200],'EdgeColor','k');
%         [C,h] = contour(XX_h/1000,YY_h/1000,hhb,[-3500:500:-500],'EdgeColor','k');
        clabel(C,h);
        hold off;
        colorbar;
        colormap(cmocean('balance'));
        caxis([-0.4 0.4]);        
        title(strcat(['$Ro=\zeta/|f0|$ at $t=$',num2str(t/t1day,'%.2f'),' days']),'interpreter','latex');        
        xlabel('$x$ (km)','interpreter','latex');
        ylabel('$y$ (km)','interpreter','latex');                   
        set(gca,'FontSize',14);
        
        
      %%% Plot zonally-averaged meridional velocity
      case 'va'
        
        %%% Load v
        data_file = fullfile(dirpath,[OUTN_V,num2str(layer-1),'_n=',num2str(n),'.dat']);
        vv = readOutputFile(data_file,Nx,Ny);    
            
        %%% Make the plot
        plot(yy_v/1000,mean(vv,1));        
        hold on
        plot(yy_v/1000,0*yy_v,'k--');
        hold off
        title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']));        
        axis([0 Ly/1000 -.02 .02]);
    
      %%% Zonally-averaged isopycnals
      case 'i'

        %%% Calculate layer surface height and load zonal velocity        
        eta = zeros(Nx,Ny,Nlay+1);
        eta(:,:,Nlay+1) = hhb;
        uu = zeros(Nx,Ny,Nlay);
        for k=Nlay:-1:1

          %%% Load kth layer thickness
          data_file = fullfile(dirpath,[OUTN_H,num2str(k-1),'_n=',num2str(n),'.dat']);
          hh = readOutputFile(data_file,Nx,Ny);    
          
          data_file = fullfile(dirpath,[OUTN_U,num2str(k-1),'_n=',num2str(n),'.dat']);
          uu(:,:,k) = readOutputFile(data_file,Nx,Ny);     
          
          %%% Add layer thickness to eta
          eta(:,:,k) = eta(:,:,k+1) + hh;
          
        end                
        
        %%% Create grid for slice contour plots
        eps = 1; %%% Small distance by which to perturb plotting points away from actual isopycnal elevations
        slice_idx = Nx/2;
        % slice_idx = 3*Nx/8;
        ZZ_slice = zeros(Ny,2*Nlay);
        UU_slice = zeros(Ny,2*Nlay);
        for k=1:Nlay
          ZZ_slice(:,2*k-1) = squeeze(eta(slice_idx,:,k))-eps;
          ZZ_slice(:,2*k) = squeeze(eta(slice_idx,:,k+1))+eps;
          UU_slice(:,2*k-1) = squeeze(uu(slice_idx,:,k));
          UU_slice(:,2*k) = squeeze(uu(slice_idx,:,k));
        end
        YY_slice = repmat(yy_h',[1 2*Nlay]);
        
        %%% Make the plot
        pcolor(YY_slice,ZZ_slice,UU_slice);
        shading interp
        hold on;
        plot(yy_h,squeeze(eta(slice_idx,:,:)),'Color',[.7 .7 .7]);
        hold off;
        colormap redblue(40)
        colorbar;
        caxis([-.2 .2])
        set(gca,'XDir','reverse');
    
%         %%% Make the plot
% %         plot(yy_h, squeeze(mean(eta,1)));     
% %         plot(yy_h,squeeze(eta(Nx/2,:,:)));
%         plot(yy_h,squeeze(eta(Nx/4,:,:)));
%         title(strcat(['t=',num2str(t/t1day,'%.2f'),' days']));        
%         xlabel('x');
%         ylabel('z');        
        
      %%% 3D isopycnals
      case 'i3d'

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
          data_file = fullfile(dirpath,[OUTN_U,num2str(layer-1),'_n=',num2str(n),'.dat']);
          uu = readOutputFile(data_file,Nx,Ny);
          data_file = fullfile(dirpath,[OUTN_V,num2str(layer-1),'_n=',num2str(n),'.dat']);
          vv = readOutputFile(data_file,Nx,Ny);           
          zeta(2:Nx,2:Ny,k) = (uu(1:Nx-1,1:Ny-1)-uu(1:Nx-1,2:Ny)) / dy + (vv(2:Nx,2:Ny)-vv(1:Nx-1,2:Ny)) / dx;                      
          
        end

        %%% Correct boundaries
        if (~useWallNS)
          zeta(:,Ny+1,:) = zeta(:,1,:);
        end
        zeta(Nx+1,:,:) = zeta(1,:,:);        
            
        
        %%% Make the plot
        clf;
        pbaspect([2,1,1]);
        eta_plot = eta(:,:,2);
        Ro = zeta(:,:,2)./(2*Omega_z);
        p = surface(XX_h/1000,YY_h/1000,eta_plot,Ro);
        p.FaceColor = 'texturemap';
        colormap(cmocean('balance'));
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
        title(strcat(['t=',num2str(t/t1year,'%.2f'),' years']));        
        xlabel('x (km)');
        ylabel('y (km)');
        zlabel('z (m)');
        set(gca,'FontSize',16);
        set(gcf,'Color','w');        
        set(gca,'ZLim',[-4650 0]);
    end
    
             
    %%% Store the image in the movie buffer
    nextframe = getframe(gcf);        
    M(Mcnt) = nextframe;       
    Mcnt = Mcnt + 1;
    
  end
  
end
