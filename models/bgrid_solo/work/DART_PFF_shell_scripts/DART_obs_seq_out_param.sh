#!/bin/bash
#
# the parameters for modifying the obs_seq.out files

# Need to look at an example obs_seq.out file in order to modify
# below parameters:
export n_obs=385      # number of scalar obs
export n_time=60       # number of obs times
export obs_freq=21600 # the obs frequency [in sec]

export n_header=10   # how many lines in the header for the obs_seq.out file
export len_block=12  # how long for each obs block
export which_row=11  # which row in an obs block has the time information
export which_row1=1   # which row in an obs block has the information ( OBS XXX)
export which_row2=5   # which row in an obs block has the information ( -1 2 -1 )

export current_dir="/home/chihchi/scratch/Bgrid_project_NEW/DART_EAKF/models/bgrid_solo/work/observation"




