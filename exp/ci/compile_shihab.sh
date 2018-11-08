#!/bin/bash 
module load libs-extra
module load intel/16
module load gcc/5.3.0
module load cmake/3.7.2
module load hwloc/1.11.6-gcc-5.3.0
module load mpi-openmpi/2.1.0-gcc-5.3.0
module load starpu/1.2.1-gcc-5.3.0-openmpi

set -x 
# Check if we are already in hicma repo dir or not.
if git -C $PWD remote -v | grep -q 'https://github.com/ecrc/hicma' 
then
	# we are, lets go to the top dir (where .git is)
	until test -d $PWD/.git ;
	do
		cd ..
	done;
else
	#we are not, we need to clone the repo
	git clone https://github.com/ecrc/hicma.git
	cd hicma
fi

# Update submodules
HICMADEVDIR=$PWD 
git submodule update --init --recursive

# STARS-H
cd stars-h
mkdir -p build/installdir
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/installdir
make -j
make install
export PKG_CONFIG_PATH=$PWD/installdir/lib/pkgconfig:$PKG_CONFIG_PATH

# CHAMELEON
cd $HICMADEVDIR
cd chameleon
mkdir -p build/installdir
cd build
cmake .. -DCMAKE_BUILD_TYPE=Debug -DCHAMELEON_USE_MPI=ON  -DCMAKE_INSTALL_PREFIX=$PWD/installdir
make -j
make install
export PKG_CONFIG_PATH=$PWD/installdir/lib/pkgconfig:$PKG_CONFIG_PATH

# HICMA
cd $HICMADEVDIR
mkdir -p build/installdir
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/installdir -DHICMA_USE_MPI=ON
make -j
make install
export PKG_CONFIG_PATH=$PWD/installdir/lib/pkgconfig:$PKG_CONFIG_PATH

# TEST
cd $HICMADEVDIR
./exp/ci/test01.sh 2 - 4 4
