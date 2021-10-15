%%%
%%% create_run.m
%%%
%%% Creates a single simulation with specified control parameters.
%%%

%%% Location of runs on file system
local_home_dir = '/Volumes/Kilchoman/UCLA/Projects/AWSIM/runs';

%%% Select model configuration and parameters
config = 'rand';
grid_size = 256; %%% Default 128
wind_stress = 0.05; %%% Default 0.05
rand_force = 0.75; %%% Default 0.75
num_canyons = 4; %%% Default 4
amp_canyons = 25; %%% Default 25
max_slope = 0.15; %%% Default 0.15
sb_width = 5; %%% Default 5
baro_force = 0; %%% Default 0
drag_coeff = 2; %%% Default 2

%%% Generate simulation name
run_name = constructRunName (config,grid_size,wind_stress, ...
          rand_force,num_canyons,amp_canyons,max_slope,sb_width,baro_force,drag_coeff);
        
%%% Create simulation
setparams (local_home_dir,run_name,config,grid_size,wind_stress, ...
          rand_force,num_canyons,amp_canyons,max_slope,sb_width,baro_force,drag_coeff);