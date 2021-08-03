%%%
%%% calc_iso_int
%%%
%%% Convenience function to calculate along-isobath integrals on isobaths specified by the vector dd. 
%%% x and y gridspacings are dx and dy. 
%%% Assumes that the curl (Nx x Ny x Nlay) and the bathymetric elevation
%%% etab (Nx x Ny) are on vorticity points.
%%%
function data_int = calc_iso_int (data,dx,dy,dd,etab)

  %%% Grid refinement factor
  ffac = 4;
  dx_f = dx/ffac;
  dy_f = dy/ffac;
  
  %%% Grid dimensions  
  Nx = size(data,1);
  Ny = size(data,2);
  Nlay = size(data,3);
  Nd = length(dd);
  Lx = Nx*dx;
  Ly = Ny*dy;  
  
  %%% Grids for interpolation
  xx = 0:dx:Lx;
  yy = 0:dy:Ly-dy;
  xx_f = 0:dx_f:Lx-dx_f; %%% N.B. don't inclue point at x=Lx to avoid double counting on periodic boundary condition
  yy_f = 0:dy_f:Ly-dy; %%% N.B. we ignore y in [Ly-dy,Ly] because it doesn't contribute to along-slope flows
  Nx_f = length(xx_f);
  Ny_f = length(yy_f);
  [XX,YY] = meshgrid(xx,yy);
  [XX_f,YY_f] = meshgrid(xx_f,yy_f);
  
  %%% Extend input grid to account for periodicity
  data = [data ; data(1,:,:)];  
  etab = [etab ; etab(1,:,:)];
  
  %%% Interpolate data onto fine grid
  data_f = zeros(Nx_f,Ny_f,Nlay);
  for k=1:Nlay
    data_f(:,:,k) = interp2(XX,YY,data(:,:,k)',XX_f,YY_f,'linear')';
  end
  etab_f = interp2(XX,YY,etab',XX_f,YY_f,'linear')';

  %%% Calculate contour lengths
  data_int = zeros(Nd,Nlay);
  for n=1:Nd
    msk_dxc = repmat( xor(-etab_f(1:Nx_f,1:Ny_f-1)>dd(n),-etab_f(1:Nx_f,2:Ny_f)>dd(n)) , [1 1 Nlay]);
    msk_dyc = repmat( xor(-etab_f(1:Nx_f,1:Ny_f)>dd(n),-etab_f([2:Nx_f 1],1:Ny_f)>dd(n)), [1 1 Nlay]);        
    data_dxc = 0.5*(data_f(:,1:Ny_f-1,:)+data_f(:,2:Ny_f,:));
    data_dyc = 0.5*(data_f(:,1:Ny_f,:)+data_f([2:Nx_f 1],1:Ny_f,:));    
    data_int(n,:) = sum(sum(dx_f.*msk_dxc.*data_dxc,1),2) + sum(sum(dy_f.*msk_dyc.*data_dyc,1),2);    
  end

end