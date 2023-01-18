#!/bin/bash
#
# the script creates many different obs seq files and each of them only contains one obs

# Need to look at an example obs_seq.out file in order to modify
# below parameters:
#export n_obs=3      # number of scalar obs
#export n_header=10   # how many lines in the header for the obs_seq.out file
#export len_block=12  # how long for each obs block
#export which_row=11  # which row in an obs block has the time information

#export current_dir="/home/chihchi/scratch/Bgrid_project_NEW/DART_EAKF/models/lorenz_96/work/observation"
#export sametime_file="${current_dir}/obs_seq.out_0_0"
#export modified_file="${current_dir}/obs_seq.out_0_0_seq"
export sametime_file=$1
export modified_file=$2


export n=0 # the counter

test -e obs_seq.out_modified && rm -f obs_seq.out_modified
cp $sametime_file $modified_file

for i in $(seq 0 $(($n_obs-1)) )
do

   line=`echo | awk "{print $n_header+$len_block*$i+$which_row}"`

   tt=$((10*$i))

   ss=`printf "%-6i"  $(($tt%86400))`
   dd=`printf "%-11i" $(($tt/86400))`

   echo "tt="${tt}  "ss="${ss} "dd="${dd}

   cat > script.sed << EOF
${line}c ${ss} ${dd} 
EOF

   sed -f script.sed -i $modified_file

done


