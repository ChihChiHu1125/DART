#!/bin/bash
#
# the parameters for modifying the obs_seq.out files

# Need to look at an example obs_seq.out file in order to modify
# below parameters:
export n_start=80
export n_time=1436       # number of obs times
export obs_freq=4 # the obs frequency [in sec]

tmp_prior=""
tmp_post=""
for i in $(seq $n_start $obs_freq $(($n_time)) )
do 

   tt=$((21600*$i))

   ss=`printf   $(($tt%86400))`
   dd=`printf   $(($tt/86400))`

   fn_prior="filter_prior_${dd}_${ss}.nc"
   fn_post="filter_post_${dd}_${ss}.nc"

   tmp_prior="$tmp_prior $fn_prior"
   tmp_post="$tmp_post $fn_post"

done

exp_name="expslp5_1day_1yr"
da_config="PFFli_ker_08cap_infR3eakffg_ensavgHT"

ncecat -v ps,t,u,v $tmp_prior "preassim_${exp_name}_${da_config}.nc"
ncecat -v ps,t,u,v $tmp_post "analysis_${exp_name}_${da_config}.nc"
