#!/bin/bash

#Load modules & global variables
module load payu
globalpath=`pwd`
count=0

# input path
input_path='/home/156/jm5970/short/mitgcm/input/global_particle_release'

# Loop for every initialization of the particle release:
for tt in `seq 0 308`
do
  # Create folder for running experiment.
  folder="30d_LADV_part_release_$(printf %05d ${tt%})"
  mkdir $folder
  # Modify corresponding files to setup the experiment.
  sed s-{0}-$(printf %05d ${tt%})-g config_sed.yaml > "$folder/config.yaml"
  cp ./input/* $folder/.
  sed s-input_off-'.'-g input/data.off > "$folder/data.off" 
  sed s-{0}-$(printf %05d ${tt%})-g config_sed.yaml > "$folder/config.yaml"
  sed s-{0}-30d_slice_chunk_$(printf %05d ${tt%})-g input/clear_archive.sh > "$folder/clear_archive.sh"
  cd $folder
  
  ln -s $input_path/flt_global_hex_032deg.bin $input_path/30d/30d_slice_chunk_$(printf %05d ${tt%})

  # Run the experiment. 
  payu run -i 0  

  cd $globalpath

  count=$((count+1))
  # Sleep for 8 hours so the process can be executed, without over queuing PBS.
  if [ $count -eq 3 ]
    then
      echo $count 
      sleep 8h
  fi
  
done

