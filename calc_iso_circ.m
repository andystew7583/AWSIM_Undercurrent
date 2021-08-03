%%%
%%% calc_iso_circ
%%%
%%% Convenience function to calculate along-isobath circulation on isobaths specified by the vector dd. 
%%% x and y gridspacings are dx and dy. 
%%% Assumes that the curl (Nx x Ny x Nlay) and the bathymetric elevation
%%% etab (Nx x Ny) are on vorticity points.
%%%
function circ = calc_iso_circ (curl,dx,dy,dd,etab)

  %%% Grid refinement factor
  ffac = 4;
  dx_f = dx/ffac;
  dy_f = dy/ffac;
  
  %%% Grid dimensions  
  Nx = size(curl,1);
  Ny = size(curl,2);
  Nlay = size(curl,3);
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
  curl = [curl ; curl(1,:,:)];  
  etab = [etab ; etab(1,:,:)];
  
  %%% Interpolate curl onto fine grid
  curl_f = zeros(Nx_f,Ny_f,Nlay);
  for k=1:Nlay
    curl_f(:,:,k) = interp2(XX,YY,curl(:,:,k)',XX_f,YY_f,'linear')';
  end
  etab_f = interp2(XX,YY,etab',XX_f,YY_f,'linear')';
  
  %%% Perform area intergrals to compute circulation
  circ = zeros(Nd,Nlay);
  for n = 1:Nd
    msk = repmat(etab_f>-dd(n),[1 1 Nlay]);    
    circ(n,:) = -sum(sum(curl_f.*msk*dx_f*dy_f,1),2);    
  end
  
end