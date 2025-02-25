#!/usr/bin/env bash

# DART software - Copyright UCAR. This open source software is provided
# by UCAR, "as is", without charge, subject to all terms of use at
# http://www.image.ucar.edu/DAReS/DART/DART_download

main() {

export DART=$(git rev-parse --show-toplevel)
source "$DART"/build_templates/buildfunctions.sh

MODEL=clm
LOCATION=threed_sphere


programs=(
filter
model_mod_check
perfect_model_obs
perturb_single_instance
)

serial_programs=(
advance_time
create_fixed_network_seq
create_obs_sequence
obs_diag
obs_seq_to_netcdf
obs_sequence_tool
)

model_serial_programs=(
clm_to_dart
dart_to_clm
)

arguments "$@"

# clean the directory
\rm -f -- *.o *.mod Makefile .cppdefs

# build and run preprocess before making any other DART executables
buildpreprocess

# build DART
buildit

# clean up
\rm -f -- *.o *.mod

}

main "$@"
