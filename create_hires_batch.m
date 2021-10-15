%%%
%%% create_run.m
%%%
%%% Creates a single simulation with specified control parameters.
%%%

%%% Load constant parameters
constants;

%%% Location of runs on file system
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';

%%% Load parameters for all simulations in the batch
[is_wind_batch,grid_size_batch,wind_stress_batch, ...
  rand_force_batch,num_canyons_batch,amp_canyons_batch, ...
  max_slope_batch,sb_width_batch,baro_force_batch, ...
  drag_coeff_batch,end_time_batch] = loadBatchParams ();

%%% Loop over simulations and create hi-res runs
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
  
  %%% Ignore existing hi-res runs
  if (grid_size ~= 128)
    continue;
  end
  
  %%% Generate simulation name
  run_name = constructRunName (config,grid_size,wind_stress, ...
          rand_force,num_canyons,amp_canyons,max_slope,sb_width,baro_force,drag_coeff);
  disp(['Working on ',run_name,' ...']);
        
  %%% Create simulation
  setparams (local_home_dir,run_name,config,grid_size,wind_stress, ...
          rand_force,num_canyons,amp_canyons,max_slope,sb_width,baro_force,drag_coeff);

end