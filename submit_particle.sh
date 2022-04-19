#!/bin/bash

#Load modules & global variables
module use /g/data3/hh5/public/modules
module load conda/analysis3-unstable

globalpath=`pwd`
count=0
cc=0
n=30
nr=729
particle_grid='flt_global_hex_032deg.bin'
freq=30
expt="${freq}d"
expt="1d"
walltime=3 #hrs


# input path
input_path='/scratch/x77/jm5970/mitgcm/input/global_particle_release'

# Loop for every initialization of the particle release:
for tt in `seq 11 $nr`
do
  # Create folder for running experiment.
  folder="${expt}/${expt}_LADV_part_release_$(printf %05d ${tt%})"
  mkdir $folder
  # Modify corresponding files to setup the experiment.
  sed -e "s-{0}-$(printf %05d ${tt%})-g;s-{1}-$expt-g;s-{2}-$(printf %02d ${walltime})-g" config_sed.yaml > "$folder/config.yaml"
  cp ./input/* $folder/.
  sed -e "s-{0}-${freq}-g" input/data > "$folder/data" 
  sed -e "s-{0}-.-g;s-{1}-$(((freq-1)*86400))-g" input/data.off > "$folder/data.off" 
  sed s-flt_global_hex_10deg.bin-${particle_grid}-g input/data.flt > "$folder/data.flt" 
  sed -e "s-{0}-${expt}_slice_chunk_$(printf %05d ${tt%})-g;s-{1}-$expt-g" input/clear_archive.sh > "$folder/clear_archive.sh"
  cd $folder
  
  ln -s $input_path/${particle_grid} $input_path/${expt}/${expt}_slice_chunk_$(printf %05d ${tt%})/

  # Initiate payu setup
  payu sweep
  payu setup && payu sweep
  
  # Run the experiment. 
  payu run -i 0  

  cd $globalpath

  count=$((count+1))

  # Sleep for 8 hours so the process can be executed, without over queuing PBS.
  if [ $cc -eq $n ]
    then
      cc=0
      echo "Sleep submission"
      sleep 2h
  else
     cc=$((cc+1))
  fi
  
done

