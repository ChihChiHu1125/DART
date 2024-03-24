#!/bin/bash
#
# a shell script that generates the directory that is necessary for DART_cylcing

mkdir state_output
mkdir obs_space_diag
mkdir log_files
mkdir temp_for_input_nml

cp filter_prior_20_0.nc ./state_output
cp merge_netcdf_files.sh ./state_output
cp ./matlab_diag_scripts/state_diag_noplot.m ./state_output
cp ./matlab_diag_scripts/cycling_read_obs_seq_files.m ./obs_space_diag
cp ./matlab_diag_scripts/obs_diag_state_RMSE_spread_at_obs_loc_noplot.m ./obs_space_diag
