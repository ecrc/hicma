export LC_ALL=en_US.UTF-8
export CRAYPE_LINK_TYPE=dynamic
#module switch PrgEnv-cray/5.2.82 PrgEnv-gnu
module switch PrgEnv-cray PrgEnv-gnu
module load hwloc/1.11.9

#module load starpu/1.2.3


#module unload PrgEnv-cray/5.2.82 
#module unload PrgEnv-intel 
#module load PrgEnv-gnu
module unload intel
#module load intel/16.3.3.210
module load intel
module load gsl
module load python  #for compiling starsh

