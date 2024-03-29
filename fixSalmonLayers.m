%%%
%%% fixSalmonLayers.m
%%%
%%% Takes output from RLSWchannel and adjusts layer thicknesses so that 
%%% they are compatible with the specified Salmon layer thickness. Must be
%%% run after regridOutput.
%%%
%%% src_home_dir: Path to directory containing the source simulation
%%% src_run_name: Name of the source simulation
%%% srcIter: Output file index from which to read the source simulation state
%%% src_home_dir: Path to directory containing the destination simulation
%%% src_run_name: Name of the destination simulation
%%% destIter: Output file index to which to write the corrected layer thicknesses
%%%
function fixSalmonLayers (src_home_dir, src_run_name, srcIter, dest_home_dir, dest_run_name, destIter)

  %%% Load source simulation parameters
  local_home_dir = src_home_dir;
  run_name = src_run_name;  
  loadParams;
  h0_src = h0;
  gsum_src = cumsum(gg);
  gg_src = gg;
  XX_src = XX_h;
  YY_src = YY_h;
  Nx_src = Nx;
  Ny_src = Ny;
  Nlay_src = Nlay;
  dirpath_src = dirpath;
  etab_src = zeros(Nx+2,Ny+2);
  etab_src(2:Nx+1,2:Ny+1) = hhb;
  [xx_src yy_src XX_src YY_src] = createmesh(-0.5*dx,Lx+0.5*dx,Nx_src+2,-0.5*dy,Ly+0.5*dy,Ny_src+2); %%% Extended grid with ghost points for interpolation
  
  %%% Load source layer thicknesses and surface pressure
  hh_src = zeros(Nlay_src,Nx_src+2,Ny_src+2);  
  MM_src = zeros(Nlay_src,Nx_src+2,Ny_src+2);  
  data_file = fullfile(dirpath_src,createOutputFilename(OUTN_PI,srcIter,-1));
  MM_src(1,2:Nx+1,2:Ny+1) = readOutputFile(data_file,Nx_src,Ny_src);      
  for k=1:Nlay        
    data_file = fullfile(dirpath_src,createOutputFilename(OUTN_H,srcIter,k));
    hh_src(k,2:Nx+1,2:Ny+1) = readOutputFile(data_file,Nx_src,Ny_src);    
  end
  
   %%% Fill in ghost points depending on boundary conditions
  if (useWallNS)    
    etab_src(:,1) = etab_src(:,2);
    etab_src(:,Ny+2) = etab_src(:,Ny+1);
    MM_src(1,:,1) = MM_src(1,:,2);
    MM_src(1,:,Ny+2) = MM_src(1,:,Ny+1);
    hh_src(:,:,1) = hh_src(:,:,2);
    hh_src(:,:,Ny+2) = hh_src(:,:,Ny+1);    
  else  
    etab_src(:,1) = etab_src(:,Ny+1);
    etab_src(:,Ny+2) = etab_src(:,2);
    MM_src(1,:,1) = MM_src(1,:,Ny+1);
    MM_src(1,:,Ny+2) = MM_src(1,:,2);
    hh_src(:,:,1) = hh_src(:,:,Ny+1);
    hh_src(:,:,Ny+2) = hh_src(:,:,2);
  end    
  if (useWallEW)    
    etab_src(1,:) = etab_src(2,:);
    etab_src(Nx+2,:) = etab_src(Nx+1,:);
    MM_src(1,1,:) = MM_src(1,2,:);
    MM_src(1,Nx+2,:) = MM_src(1,Nx+1,:);
    hh_src(:,1,:) = hh_src(:,2,:);
    hh_src(:,Nx+2,:) = hh_src(:,Nx+1,:);
  else    
    etab_src(1,:) = etab_src(Nx+1,:);
    etab_src(Nx+2,:) = etab_src(2,:);
    MM_src(1,1,:) = MM_src(1,Nx+1,:);
    MM_src(1,Nx+2,:) = MM_src(1,2,:);
    hh_src(:,1,:) = hh_src(:,Nx+1,:);
    hh_src(:,Nx+2,:) = hh_src(:,2,:);
  end   
  
  %%% Construct source layer elevations
  ee_src = zeros(Nlay_src+1,Nx_src+2,Ny_src+2);
  ee_src(2:Nlay_src+1,:,:) = -cumsum(hh_src);
  ee_src(Nlay_src+1,:,:) = etab_src;
   
  %%% Construct source Montgomery potential everywhere
  for k = 2:Nlay_src
    MM_src(k,:,:) = MM_src(k-1,:,:) + gg(k)*ee_src(k,:,:) - gsum_src(k)/3 * h0_src.^4 ./ hh_src(k,:,:).^3;
  end
    
  %%% Load destination simulation parameters
  local_home_dir = dest_home_dir;
  run_name = dest_run_name;
  loadParams;
  h0_dest = h0;
  gsum_dest = cumsum(gg);
  gg_dest = gg;
  XX_dest = XX_h;
  YY_dest = YY_h;
  Nx_dest = Nx;
  Ny_dest = Ny;
  Nlay_dest = Nlay;
  dirpath_dest = dirpath;
  
  
  %%% Interpolate Montgomery potential onto destination grid
  hh_dest = zeros(Nlay_dest,Nx_dest,Ny_dest);
  MM_dest = zeros(Nlay_dest,Nx_dest,Ny_dest);
  for k=1:Nlay
    hh_dest(k,:,:) = interp2(XX_src',YY_src',squeeze(hh_src(k,:,:))',XX_dest',YY_dest','linear')';    
    MM_dest(k,:,:) = interp2(XX_src',YY_src',squeeze(MM_src(k,:,:))',XX_dest',YY_dest','linear')';    
  end
  etab_dest = interp2(XX_src',YY_src',squeeze(etab_src(:,:))',XX_dest',YY_dest','linear')';
  
  %%% Adjust layer thicknesses to account for incrops
  hh_dest = solveSalmonThicknesses (hh_dest,etab_dest,zeros(Nx_dest,Ny_dest),MM_dest,gg_dest,h0_dest); 
  
  %%% Write corrected layer thickensses to output file
  for k=1:Nlay
    fname = fullfile(dirpath_dest,createOutputFilename(OUTN_H,destIter,k));
    writeDataFile(fname,hh_dest(k,:,:));
  end
  
  %%% Write topography
  fname = fullfile(dirpath_dest,'etab.dat');
  writeDataFile(fname,etab_dest);
    
end