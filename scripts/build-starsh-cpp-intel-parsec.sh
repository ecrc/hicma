#!/bin/bash -le
set -x
currdir=$PWD
function module_load {
	module=$1
	if module list 2>&1 | grep $module; then 
		echo "$module is loaded"; 
	else
		module load $module
	fi
}
function download {
	reponame=$1
	# Check if we are already in repo dir or not.
	if git -C $PWD remote -v | grep -q "https://github.com/ecrc/$reponame"
	then
		echo "We are in $reponame"
		# we are, lets go to the top dir (where .git is)
		until test -d $PWD/.git ;
		do
			cd ..
		done;
	else
		echo "We are NOT in $reponame"
		if [ ! -d $reponame ]; then
			#we are not, we need to clone the repo
			git clone git@github.com:ecrc/$reponame.git
			if [ ! -d $reponame ]; then
				echo "Failed to clone $reponame"
				return
			fi
		fi
		cd $reponame
	fi
}


install_starsh_core=1
install_starsh=1

module_load cmake/3.11.1
module_load intel/2018

if [ $install_starsh_core -eq 1 ]; then
    reponame=stars-h-core-dev
    download $reponame
    git submodule update --init
    if [ -d build-intel ]; then
        rm -rf build-intel
    fi
    mkdir build-intel
    cd build-intel
    cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/install
    make -j install
    export LD_LIBRARY_PATH=/opt/share/intel/2018/compilers_and_libraries/linux/lib/intel64/:$LD_LIBRARY_PATH
    #/home/akbudak/stars-h-core-dev/build-intel/tests/applications/particles
    if [ -d $PWD/install ]; then
        export PKG_CONFIG_PATH=$HOME/stars-h-core-dev/build-intel/install/lib/pkgconfig/:$PKG_CONFIG_PATH
    fi
    cd $currdir
fi

if [ $install_starsh -eq 1 ]; then
    module_load hwloc/1.11.8-intel-2018
    module_load plasma/2.8.0-intel-2018-mkl 
    module_load parsec/master-intel-2018-mkl-intelmpi-plasma-2.8.0
    reponame=stars-h-dev
    download $reponame
    git submodule update --init
    git checkout muxas/cpp
    if [ -d build-intel ]; then
        rm -rf build-intel
    fi
    mkdir build-intel
    cd build-intel
    cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/install -DPARSEC_daguepp_BIN_DIR=/opt/ecrc/parsec/master-intel-2018-mkl-intelmpi-plasma-2.8.0/ub16/bin -DPARSEC_SRC_DIR=$HOME/parsec
    make -j install
    if [ -d $PWD/install ]; then
        export PKG_CONFIG_PATH=$HOME/stars-h-dev/build-intel/install/lib/pkgconfig/:$PKG_CONFIG_PATH
    fi
    export MKL_NUM_THREADS=1 #starsh is linked to parallel mkl 
    cd $currdir
fi
