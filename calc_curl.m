%%%
%%% calc_curl
%%%
%%% Convenience function to calculate curl on vorticity points.
%%%
function curl = calc_curl (u,v,dx,dy)

  Nx = size(u,1);
  Ny = size(u,2);
  Nlay = size(u,3);  
  
  %%% Compute curl, including points on southern boundary
  curl = zeros(Nx,Ny,Nlay);
  curl(:,2:Ny,:) = (v([1:Nx],2:Ny,:)-v([Nx 1:Nx-1],2:Ny,:))/dx - (u(:,2:Ny,:)-u(:,1:Ny-1,:))/dy;
  curl(:,1,:) = (v([1:Nx],1,:)-v([Nx 1:Nx-1],1,:))/dx - u(:,1,:)/dy;
  
end