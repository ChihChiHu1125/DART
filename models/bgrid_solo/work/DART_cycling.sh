#!/bin/bash
#
# a shell script that assimilates obs and then advances the ensemble in time
# Chih-Chi Hu 2022/12/30

export obs_freq=21600 # obs frequency in seconds (check the folder ./observation for info there)
export n_da_cycle=60   # how many DA cycles to run (note: 1st DA cycle starts from 0 day 0 sec)

export current_dir="/home/chihchi/scratch/Bgrid_project_NEW/DART_EAKF/models/lorenz_96/work"
export obs_files="observation"
export obs_space_diag="obs_space_diag"
export state_space_diag="state_output"
export log_file="log_files"


n_cycle=0 # the counter for DA cycle

for i in $(seq 0 $((n_da_cycle-1)) )
do
   # time information for this DA cycle:
   tt=$(($obs_freq*$i))
   ss="$(($tt%86400))"
   dd="$(($tt/86400))"

   tt_prev=$(($obs_freq*($i-1)))
   ss_prev="$(($tt_prev%86400))"
   dd_prev="$(($tt_prev/86400))"

   echo "start doing DA cycle " ${i} "time = "${dd} "day" ${ss} "second."

   if [ "$i" -eq 0 ]; then

      # STEP 1: ensemble forecasting and record the obs space prior
      echo "=== STEP 1: ensemble forecasting and record the obs space prior ==="

      cp input.nml.template input.nml

      state_input_file="${state_space_diag}\/filter_prior_${dd}_${ss}.nc"
      state_output_file="${state_space_diag}\/useless.nc"
      stage_to_write="'null'"
      obs_input_file="${obs_files}\/obs_seq.out_0_0"
      obs_output_file="${obs_space_diag}\/obs_seq.prior_0_0"
      evaluate_or_assimilate="evaluate_these_obs_types"
      async="0"                                            # if evaluate, then async = 0
      init_time_days="0"
      init_time_seconds="0"

      cat > script.sed << EOF
         s/cchu_input_file_here.nc/${state_input_file}/g;
         s/cchu_output_file_here.nc/${state_output_file}/g;
         s/stages_to_write  =.*/stages_to_write  = ${stage_to_write}/g;
         s/cchu_obs_seq.out/${obs_input_file}/g;
         s/cchu_obs_seq.final/${obs_output_file}/g;
         s/assimilate_these_obs_types/${evaluate_or_assimilate}/g;
         s/async                        =.*/async                        = ${async},/g;
         s/init_time_days               =.*/init_time_days               = ${init_time_days},/g;
         s/init_time_seconds            =.*/init_time_seconds            = ${init_time_seconds},/g
EOF
      sed -f script.sed -i input.nml

      cp input.nml ./temp_for_input_nml/"input.nml_prior_${dd}_${ss}" # save the input.nml files for debugging

      ./filter >& ${log_file}/prior_log_${dd}_${ss}
      #mpirun -n 5 ./filter >& ${log_file}/prior_log_${dd}_${ss}

    else # if not the first DA cycle:

      # STEP 1: ensemble forecasting and record the obs space prior
      echo "=== STEP 1: ensemble forecasting and record the obs space prior ==="

      cp input.nml.template input.nml

      state_input_file="${state_space_diag}\/filter_post_${dd_prev}_${ss_prev}.nc"
      state_output_file="${state_space_diag}\/filter_prior_${dd}_${ss}.nc"
      stage_to_write="'output'"
      obs_input_file="${obs_files}\/obs_seq.out_${dd}_${ss}"
      obs_output_file="${obs_space_diag}\/obs_seq.prior_${dd}_${ss}"
      evaluate_or_assimilate="evaluate_these_obs_types"
      async="0"                                            # if evaluate, then async = 0
      init_time_days="${dd_prev}"
      init_time_seconds="${ss_prev}"

      cat > script.sed << EOF
         s/cchu_input_file_here.nc/${state_input_file}/g;
         s/cchu_output_file_here.nc/${state_output_file}/g;
         s/stages_to_write  =.*/stages_to_write  = ${stage_to_write}/g;
         s/cchu_obs_seq.out/${obs_input_file}/g;
         s/cchu_obs_seq.final/${obs_output_file}/g;
         s/assimilate_these_obs_types/${evaluate_or_assimilate}/g;
         s/async                        =.*/async                        = ${async},/g;
         s/init_time_days               =.*/init_time_days               = ${init_time_days},/g;
         s/init_time_seconds            =.*/init_time_seconds            = ${init_time_seconds},/g
EOF

      sed -f script.sed -i input.nml

      cp input.nml ./temp_for_input_nml/"input.nml_prior_${dd}_${ss}"
      #./filter >& ${log_file}/prior_log_${dd}_${ss}
      mpirun -n 5 ./filter >& ${log_file}/prior_log_${dd}_${ss}

   fi


   # STEP 2: DA step
   echo "=== STEP 2: sequentially assimilating scalar obs                ==="

   cp input.nml.template input.nml

   state_input_file="${state_space_diag}\/filter_prior_${dd}_${ss}.nc"
   state_output_file="${state_space_diag}\/filter_post_${dd}_${ss}.nc"
   stage_to_write="'output'"
   obs_input_file="${obs_files}\/obs_seq.out_${dd}_${ss}_seq"
   obs_output_file="${obs_space_diag}\/useless"
   evaluate_or_assimilate="assimilate_these_obs_types"
   async="1"                                            # if assimilate, then async = 1
   init_time_days="0"                                   # if assimilate, set init time = 0 day ,0 sec
   init_time_seconds="0"

   cat > script.sed << EOF
      s/cchu_input_file_here.nc/${state_input_file}/g;
      s/cchu_output_file_here.nc/${state_output_file}/g;
      s/stages_to_write  =.*/stages_to_write  = ${stage_to_write}/g;
      s/cchu_obs_seq.out/${obs_input_file}/g;
      s/cchu_obs_seq.final/${obs_output_file}/g;
      s/assimilate_these_obs_types/${evaluate_or_assimilate}/g;
      s/async                        =.*/async                        = ${async},/g;
      s/init_time_days               =.*/init_time_days               = ${init_time_days},/g;
      s/init_time_seconds            =.*/init_time_seconds            = ${init_time_seconds},/g;
      s/compute_posterior            =.*/compute_posterior            = .false.,/g;
      s/dt_atmos =.*/dt_atmos = 10,/g
EOF
   sed -f script.sed -i input.nml

   cp input.nml ./temp_for_input_nml/"input.nml_DA_${dd}_${ss}"
   #./filter >& ${log_file}/DA_log_${dd}_${ss}
   mpirun -n 5 ./filter >& ${log_file}/DA_log_${dd}_${ss}

   # STEP 3: Record the obs space diagnostics for posterior
   echo "=== STEP 3: Record the obs space diagnostics for posterior      ==="

   cp input.nml.template input.nml

   state_input_file="${state_space_diag}\/filter_post_${dd}_${ss}.nc"
   state_output_file="${state_space_diag}\/useless.nc"
   stage_to_write="'null'"
   obs_input_file="${obs_files}\/obs_seq.out_${dd}_${ss}"
   obs_output_file="${obs_space_diag}\/obs_seq.post_${dd}_${ss}"
   evaluate_or_assimilate="evaluate_these_obs_types"
   async="0"                                            # if assimilate, then async = 1
   init_time_days="${dd}"                                   # if assimilate, set init time = 0 day ,0 sec
   init_time_seconds="${ss}"

   cat > script.sed << EOF
      s/cchu_input_file_here.nc/${state_input_file}/g;
      s/cchu_output_file_here.nc/${state_output_file}/g;
      s/stages_to_write  =.*/stages_to_write  = ${stage_to_write}/g;
      s/cchu_obs_seq.out/${obs_input_file}/g;
      s/cchu_obs_seq.final/${obs_output_file}/g;
      s/assimilate_these_obs_types/${evaluate_or_assimilate}/g;
      s/async                        =.*/async                        = ${async},/g;
      s/init_time_days               =.*/init_time_days               = ${init_time_days},/g;
      s/init_time_seconds            =.*/init_time_seconds            = ${init_time_seconds},/g
EOF
   sed -f script.sed -i input.nml

   cp input.nml ./temp_for_input_nml/"input.nml_post_${dd}_${ss}"
   ./filter >& ${log_file}/post_log_${dd}_${ss}
   #mpirun -n 5 ./filter >& ${log_file}/post_log_${dd}_${ss}

done

echo " === !!! All cycling experiments done !!! === "

