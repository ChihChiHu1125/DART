#!/bin/bash
#
# the script separates an obs_seq.out file which contains many obs times into 
# several files which each of them contain only one obs time.

# Need to look at an example obs_seq.out file in order to modify
# below parameters:
#export n_obs=3        # number of scalar obs
#export n_time=1       # number of obs times
#export obs_freq=43200 # the obs frequency [in sec]
#export n_header=10    # how many lines in the header for the obs_seq.out file
#export len_block=12   # how long for each obs block
#export which_row1=1   # which row in an obs block has the information that needs to be changed
#export which_row2=5   # which row in an obs block has the information that needs to be changed
#export current_dir="/home/chihchi/scratch/Bgrid_project_NEW/DART_EAKF/models/lorenz_96/work/observation"
export all_time_file=$1

source DART_obs_seq_out_param.sh

n=0 # the counter

while [ $n -lt $n_time ];
do

   tt=$(($obs_freq*$n))

   ss=`printf $(($tt%86400))`
   dd=`printf $(($tt/86400))`

   echo "the ${n}-th obs"  "ss="${ss} "dd="${dd}

   modified_file=${1}_${dd}_${ss}

   test -e $modified_file && rm -f $modified_file
   cp $all_time_file $modified_file

   # STEP 1: modify the header info:

   # replace the header info, and delete the actual obs info:
   cat > script.sed << EOF
      6c  \  num_obs:            $n_obs  max_num_obs:            $n_obs
      10c \  first:            1  last:            $n_obs
      11,\$d
EOF

   sed -f script.sed -i $modified_file

   # STEP 2: add back the obs info at the required time:
   needed_line_start=$(( $n_header + $n*$len_block*n_obs + 1 ))
   needed_line_end=$((   $n_header + ($n+1)*$len_block*n_obs ))

   cat > script.sed << EOF
      ${needed_line_start},${needed_line_end}p
EOF

   sed -f script.sed -n $all_time_file >> $modified_file

   # STEP 3: modify each obs info:
   for i in $(seq 1 $n_obs )
   do
      modify_line1=$(( $n_header + $(($i-1))*$len_block + $which_row1 ))
      modify_line2=$(( $n_header + $(($i-1))*$len_block + $which_row2 ))

      if [ $i == 1 ]; then
         cat > script.sed << EOF
            ${modify_line1}c  \ OBS            $(($i))
            ${modify_line2}c  \          -1           2          -1
EOF

         sed -f script.sed -i $modified_file

      elif [ $i == $n_obs ]; then
         cat > script.sed << EOF
            ${modify_line1}c  \ OBS            $(($i))
            ${modify_line2}c  \           $(($i-1))          -1          -1
EOF

         sed -f script.sed -i $modified_file

      else # i != 1 or n_obs
         cat > script.sed << EOF
            ${modify_line1}c  \ OBS            $(($i))
            ${modify_line2}c  \           $(($i-1))           $(($i+1))          -1
EOF

         sed -f script.sed -i $modified_file

      fi
   done

   seq_file=${modified_file}_seq

   ./DART_modify_separate_obs_seq_out.sh $modified_file ${seq_file}

   n=$(($n+1))

done # while loop
