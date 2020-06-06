#conda deactivate
module load cmake/3.13.4
export LC_ALL=en_US.UTF-8
module load intel/19.0.5.281
module load cray-mpich/7.7.11
module load python
module switch PrgEnv-cray PrgEnv-gnu
module unload cray-libsci
module list -l
export CRAYPE_LINK_TYPE=dynamic


export PKG_CONFIG_PATH=/project/k1205/omairyrm/starpu-1.2.6/install/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=/project/k1205/omairyrm/starpu-1.2.6/install/lib:$LD_LIBRARY_PATH
#hwloc
export HWLOC_SRC_DIR=/project/k1205/omairyrm/hwloc/hwloc-1.11.10
export HWLOC_ROOT=${HWLOC_SRC_DIR}/install
export PKG_CONFIG_PATH=${HWLOC_ROOT}/lib/pkgconfig:$PKG_CONFIG_PATH
export LD_LIBRARY_PATH=${HWLOC_ROOT}/lib:$LD_LIBRARY_PATH


