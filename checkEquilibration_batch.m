%%%
%%% calcMomBudget_batch.m
%%%
%%% Calculates along-isobath momentum budget diagnostics for all simulations.
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

%%% Iterate through simulations 
Nsims = length(is_wind_batch);
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
        
  %%% Load time series diagnostics
  [KE,PE,E,Z,t]=readEZfile(local_home_dir,run_name);
  if (t(end)~=end_time*t1year)
    disp(['Mismatched end time for run: ',run_name]);
  end
  start_time = end_time*t1year-5*t1year;
  start_idx = find(t<=start_time,1,'last');
  
  %%% Construct a range of indices that creates a continuous time series
  idx = start_idx:length(t);  %%% Default to all indices within the time range
  m = 1;
  %%% Loop through indices within the time range
  while (m < length(idx))
    m = m + 1;
    if (t(idx(m)) < t(idx(m-1))) %%% Time goes backwards - indicates a restart      
      restart_idx = find(t(1:idx(m-1)) < t(idx(m)),1,'last'); %%% Find the previous index in t that matches the restart time
      idx_idx = find(idx == restart_idx); %%% Index within the index array of the restart time
      idx(idx_idx+1:m-1) = []; %%% Remove the chunk of time that was re-run since the restart
      m = idx_idx + 1; %%% Resume the loop, noting that now idx has changed size
    end
  end
  
  %%% Compute trends
  linfit = polyfit(t(idx),PE(idx)/mean(PE(idx)-PE(1))*100,1);
  disp(['PE trend = ',num2str(linfit(1)*5*t1year,'%.1f')]);
  linfit = polyfit(t(idx),KE(idx)/mean(KE(idx))*100,1);
  disp(['KE trend = ',num2str(linfit(1)*5*t1year,'%.1f')]);
  [r,p] = corr(t(idx)',PE(idx)');
  disp(['PE trend r-value = ',num2str(r)]);
  [r,p] = corr(t(idx)',KE(idx)');
  disp(['KE trend r-value = ',num2str(r)]);
        
end