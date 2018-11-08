#!/bin/bash -le
module load libs-extra
module load mkl/2018-initial
module load gcc/5.3.0
module load cmake/3.7.2
module load hwloc/1.11.6-gcc-5.3.0
#module load mpi-openmpi/2.1.0-gcc-5.3.0
#module load starpu/1.2.1-gcc-5.3.0-openmpi
module load starpu/1.2.1-gcc-5.3.0

module list
git config --global credential.helper 'cache --timeout=36000'

# BASH verbose mode
set -x 


reponame=hicma
# Check if we are already in hicma repo dir or not.
if git -C $PWD remote -v | grep -q "https://github.com/ecrc/$reponame"
then
	# we are, lets go to the top dir (where .git is)
	until test -d $PWD/.git ;
	do
		cd ..
	done;
else
	#we are not, we need to clone the repo
	git clone https://github.com/ecrc/$reponame.git
	cd $reponame
fi

# Update submodules
HICMADEVDIR=$PWD 
git submodule update --init --recursive

# STARS-H
cd stars-h
rm -rf build
mkdir -p build/installdir
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/installdir -DMPI=OFF -DOPENMP=OFF -DSTARPU=OFF -DBLAS_LIBDIR=$MKLROOT/lib/ -DBLAS_DIR=$MKLROOT 
make clean
make -j
make install
export PKG_CONFIG_PATH=$PWD/installdir/lib/pkgconfig:$PKG_CONFIG_PATH

# CHAMELEON
cd $HICMADEVDIR
cd chameleon
rm -rf build
mkdir -p build/installdir
cd build
cmake ..  -DCMAKE_INSTALL_PREFIX=$PWD/installdir -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=ON -DCHAMELEON_ENABLE_EXAMPLE=ON -DCHAMELEON_ENABLE_TESTING=ON -DCHAMELEON_ENABLE_TIMING=ON -DCHAMELEON_USE_MPI=OFF -DCHAMELEON_USE_CUDA=OFF -DCHAMELEON_SCHED_QUARK=OFF -DCHAMELEON_SCHED_STARPU=ON  -DSTARPU_DIR=/Users/akbudak/starpu-svn-install/ -DBLAS_LIBRARIES="-L${MKLROOT}/lib;-lmkl_intel_lp64;-lmkl_core;-lmkl_sequential;-lpthread;-lm;-ldl"
-DBLAS_COMPILER_FLAGS="-m64;-I${MKLROOT}/include" -DLAPACK_LIBRARIES="-L${MKLROOT}/lib;-lmkl_intel_lp64;-lmkl_core;-lmkl_sequential;-lpthread;-lm;-ldl" -DCBLAS_DIR="${MKLROOT}" -DLAPACKE_DIR="${MKLROOT}" -DTMG_DIR="${MKLROOT}"
make clean
make -j
make install
export PKG_CONFIG_PATH=$PWD/installdir/lib/pkgconfig:$PKG_CONFIG_PATH

# HICMA
cd $HICMADEVDIR
rm -rf build
mkdir -p build/installdir
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/installdir -DHICMA_USE_MPI=OFF -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=ON -DCHAMELEON_ENABLE_EXAMPLE=ON -DCHAMELEON_ENABLE_TESTING=ON -DCHAMELEON_ENABLE_TIMING=ON -DCHAMELEON_USE_MPI=OFF -DCHAMELEON_USE_CUDA=OFF -DCHAMELEON_SCHED_QUARK=OFF -DCHAMELEON_SCHED_STARPU=ON  -DSTARPU_DIR=/Users/akbudak/starpu-svn-install/ -DBLAS_LIBRARIES="-L${MKLROOT}/lib;-lmkl_intel_lp64;-lmkl_core;-lmkl_sequential;-lpthread;-lm;-ldl" -DBLAS_COMPILER_FLAGS="-m64;-I${MKLROOT}/include" -DLAPACK_LIBRARIES="-L${MKLROOT}/lib;-lmkl_intel_lp64;-lmkl_core;-lmkl_sequential;-lpthread;-lm;-ldl" -DCBLAS_DIR="${MKLROOT}" -DLAPACKE_DIR="${MKLROOT}" -DTMG_DIR="${MKLROOT}"
make clean
make -j
make install
export PKG_CONFIG_PATH=$PWD/installdir/lib/pkgconfig:$PKG_CONFIG_PATH

