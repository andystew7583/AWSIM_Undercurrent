%%%
%%% create_run.m
%%%
%%% Creates a single simulation with specified control parameters.
%%%

%%% Load constant parameters
constants;

%%% Location of runs on file system
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';

%%% Grid size for high-res runs vs lo-res runs
hires_grid_size = 256;
lores_grid_size = 128;
lores_iters_per_year = 73;

%%% Load parameters for all simulations in the batch
[is_wind_batch,grid_size_batch,wind_stress_batch, ...
  rand_force_batch,num_canyons_batch,amp_canyons_batch, ...
  max_slope_batch,sb_width_batch,baro_force_batch, ...
  drag_coeff_batch,end_time_batch] = loadBatchParams ();

%%% Script files
run_batch_fname = 'run_batch.sh';
run_batch_file = fopen(fullfile(local_home_dir,run_batch_fname),'w');

%%% Loop over simulations and create hi-res runs
for n=1:Nsims    
% for n=[17 19]
  
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
  
  %%% Ignore existing hi-res runs
  if (grid_size ~= lores_grid_size)
    continue;
  end
  
  %%% Generate simulation name
  hires_run_name = constructRunName (config,hires_grid_size,wind_stress, ...
          rand_force,num_canyons,amp_canyons,max_slope,sb_width,baro_force,drag_coeff);
  disp(['Working on ',hires_run_name,' ...']);
        
  %%% Create simulation
  setparams (local_home_dir,hires_run_name,config,hires_grid_size,wind_stress, ...
          rand_force,num_canyons,amp_canyons,max_slope,sb_width,baro_force,drag_coeff);
        
  %%% Name and source iteration of low-resolution run to restart from
  lores_run_name = constructRunName (config,lores_grid_size,wind_stress, ...
          rand_force,num_canyons,amp_canyons,max_slope,sb_width,baro_force,drag_coeff);        
  lores_iter = end_time*lores_iters_per_year;
  
  %%% Regrid low-resolution model state to high-resolution grid
  regridOutput(local_home_dir,lores_run_name,lores_iter,2*hires_grid_size,hires_grid_size,fullfile(local_home_dir,hires_run_name),0)
  
  %%% Regrid Salmon layers
  fixSalmonLayers(local_home_dir,lores_run_name,lores_iter,local_home_dir,hires_run_name,0)
  
  %%% Add lines to run_batch file to execute this simulation  
  fprintf(run_batch_file,'cd %s\n',hires_run_name);
  fprintf(run_batch_file,'sh Make_fftw.sh\n');
  fprintf(run_batch_file,'sh Run.sh\n');
  fprintf(run_batch_file,'cd ..\n');

end

%%% Close open file handles
fclose(run_batch_file);