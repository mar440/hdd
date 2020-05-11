#!/bin/bash

source /etc/profile.d/modules.sh
module load mpi/openmpi-x86_64
source srclibs

echo "number of args:  $#"
if (( $# < 1 )); then
    echo "add number of mpi as cmd argument"
else
  orig_path=$(pwd)

  source srclibs

  #recompile python wrapper
  rm HddApi.cpp
  python setup.py build_ext --inplace

  # launch mpi
  echo "mpisize: $1"
  mpirun -n $1 python FEM_launch.py
#  cd ${orig_path}
fi
