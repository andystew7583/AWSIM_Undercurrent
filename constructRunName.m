%%%
%%% constructRunName.m
%%%
%%% Convenience function to construct simulation names consistently.
%%%
function run_name = constructRunName (config,grid_size,wind_stress, ...
          rand_force,num_canyons,amp_canyons,max_slope,sb_width,baro_force,drag_coeff)

  %%% Load definitions
  constants;
                      
  %%% Base simulation name
  run_name = 'undercurrent';
     
  %%% Model configuration type
  switch (config)
    case 'wind'
      run_name = [run_name,'_LWF'];
      run_name = [run_name,'_tau',num2str(wind_stress)];
    case 'rand'
      run_name = [run_name,'_RF'];
      run_name = [run_name,'_F',num2str(rand_force)];
  end
  
  %%% Parameter values
  run_name = [run_name,'_Nc',num2str(num_canyons)];
  run_name = [run_name,'_Yc',num2str(amp_canyons)];
  if (max_slope ~= 0.15)
    run_name = [run_name,'_s',num2str(max_slope)];
  end
  if (sb_width ~= 5)
    run_name = [run_name,'_Ws',num2str(sb_width)];
  end
  if (baro_force ~= 0)
    run_name = [run_name,'_Fbaro',num2str(baro_force)];
  end  
  if (drag_coeff ~= 2)
    run_name = [run_name,'_Cd',num2str(drag_coeff)];
  end

  %%% Resolution
  if (grid_size ~= 128)
    run_name = [run_name,'_Ny',num2str(grid_size)];
  end
  
end