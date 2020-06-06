module swap PrgEnv-cray/6.0.5 PrgEnv-allinea
module load cdt/19.08
#export ARMPL_DIR=/opt/allinea/19.0.0/opt/arm/armpl-19.0.0_ThunderX2CN99_SUSE-12_arm-hpc-compiler_19.0_aarch64-linux
#export LD_LIBRARY_PATH=/opt/allinea/19.0.0/opt/arm/armpl-19.0.0_ThunderX2CN99_SUSE-12_arm-hpc-compiler_19.0_aarch64-linux/lib:$LD_LIBRARY_PATH
#export PKG_CONFIG_PATH=/lustre/home/ri-kakbudak/hicma-dev/stars-h/build/installdir/lib/pkgconfig:$PKG_CONFIG_PATH
ARMPLLIBS="-larmpl"
ARMPLROOT="/opt/allinea/19.2.0.0/opt/arm/armpl-19.2.0_ThunderX2CN99_SUSE-12_arm-hpc-compiler_19.2_aarch64-linux"
export LD_LIBRARY_PATH="/opt/allinea/19.2.0.0/opt/arm/armpl-19.2.0_ThunderX2CN99_SUSE-12_arm-hpc-compiler_19.2_aarch64-linux/lib":$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=/home/ri-ralomairy/sourcefiles/lapack-3.9.0/build/installdir:$LD_LIBRARY_PATH
export PKG_CONFIG_PATH=/home/ri-ralomairy/sourcefiles/lapack-3.9.0/build/installdir/pkgconfig:$PKG_CONFIG_PATH
