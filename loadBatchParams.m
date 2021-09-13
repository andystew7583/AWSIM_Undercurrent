%%%
%%% loadBatchParams.m
%%%
%%% Retrieves a tabulation of the control parameters used in all
%%% simulations.
%%%
function [is_wind_batch,grid_size_batch,wind_stress_batch, ...
          rand_force_batch,num_canyons_batch,amp_canyons_batch, ...
          max_slope_batch,sb_width_batch,baro_force_batch, ...
          drag_coeff_batch,end_time_batch] = loadBatchParams ()        

  %%% Load parameters from file
  load batchData.txt;

  %%% Extract individual parameters
  is_wind_batch = batchData(:,1);
  grid_size_batch = batchData(:,2);
  wind_stress_batch = batchData(:,3);
  rand_force_batch = batchData(:,4);
  num_canyons_batch = batchData(:,5);
  amp_canyons_batch = batchData(:,6);
  max_slope_batch = batchData(:,7);
  sb_width_batch = batchData(:,8);
  baro_force_batch = batchData(:,9);
  drag_coeff_batch = batchData(:,10);
  end_time_batch = batchData(:,11);

end