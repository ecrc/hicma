#!/bin/bash -le

MPIVALUE="OFF"
CUDAVALUE="OFF"
WORKSPACE=$HOME/hicma-dev
currdir=$PWD
build_hwloc=0
build_starpu=0
build_chameleon=0
build_starsh=0
build_starsh_core=0
python="-DPYTHON_EXECUTABLE=/usr/bin/python"
if [ ! -f /usr/bin/python ]; then
    echo "Python 2.7 is required"
    return
fi
if [ -z $MKLROOT ]; then
	echo "\$MKLROOT must be set"
	return
else 
	echo "Intel MKL library is found: $MKLROOT"
fi

git config --global credential.helper 'cache --timeout=36000'

# BASH verbose mode
set -x 


reponame=hicma-dev
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

PREFIX=${PREFIX:-$WORKSPACE/install-prefix}
TMPDIR=/tmp/.hicma-temp
BASEDIR=$WORKSPACE

export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:$PREFIX/lib/pkgconfig
export MKLROOT=${MKLROOT:-/opt/intel/mkl}
source ${MKLROOT}/mklvars.sh intel64
export DYLD_LIBRARY_PATH=$PREFIX/lib:$DYLD_LIBRARY_PATH
export CPATH=$PREFIX/include:$CPATH

#set -xe e means exit if there is error
set -x
rm -rf $TMPDIR
mkdir -p $TMPDIR
#rm -rf $PREFIX
mkdir -p $PREFIX

# lower priority of job 
renice -n 15 -p $$


# hwloc
if [ $build_hwloc -eq 1 ]; then
    echo "====================================HWLOC"
    if pkg-config --exists hwloc
    then
        _LOCATION=`pkg-config --variable=prefix hwloc`
        echo "hwloc FOUND in [$_LOCATION]"
    else
        echo "Building Hwloc..."
        cd $TMPDIR
        wget https://www.open-mpi.org/software/hwloc/v1.11/downloads/hwloc-1.11.5.tar.gz -O - | tar -zx
        cd hwloc-1.11.5
        ./configure --prefix=$PREFIX
        make -j 4 -l 5 || make && make install
    fi
fi


# StarPU
if [ $build_starpu -eq 1 ]; then
    echo "====================================StarPU"
    if pkg-config --exists --atleast-version=1.2 libstarpu
    then
        _LOCATION=`pkg-config --variable=prefix libstarpu`
        echo "StarPU FOUND in [$_LOCATION]"
    else
        echo "Building StarPU..."
        cd $TMPDIR
        wget http://starpu.gforge.inria.fr/files/starpu-1.2.5/starpu-1.2.5.tar.gz -O - | tar -zx
        cd starpu-1.2.5
        if [ "$CUDAVALUE" == "ON" ]; then
            ./configure --enable-cuda --disable-opencl --prefix=$PREFIX
        else
            #LDFLAGS="-lmpi" 
            ./configure --disable-cuda --disable-opencl --prefix=$PREFIX
        fi
        make -j 4 -l 5|| make && make install
    fi
fi

# CHAMELEON
if [ $build_chameleon -eq 1 ]; then
    echo "Building CHAMELEON..."
    cd $BASEDIR/chameleon
    rm -rf build
    mkdir -p build && cd build
    cmake .. -DCHAMELEON_USE_MPI=$MPIVALUE -DCHAMELEON_USE_CUDA=$CUDAVALUE -DCHAMELEON_ENABLE_EXAMPLE=OFF -DCHAMELEON_ENABLE_TESTING=OFF -DCHAMELEON_ENABLE_TIMING=OFF -DBUILD_SHARED_LIBS=ON -DCMAKE_INSTALL_PREFIX=$PREFIX -DBLAS_LIBRARIES="-L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl" $python 
    make -j 4 -l 5 || make VERBOSE=1
    make install
fi

# STARS-H
if [ $build_starsh -eq 1 ]; then
    echo "Building STARS-H..."
    cd $BASEDIR/stars-h
    mkdir -p build && cd build
    cmake .. -DCMAKE_C_FLAGS=-fPIC -DEXAMPLES=OFF -DTESTING=OFF -DMPI=$MPIVALUE -DCMAKE_INSTALL_PREFIX=$PREFIX
    make -j 4 -l 5|| make VERBOSE=1
    make install
fi

# STARS-H-CORE
if [ $build_starsh_core -eq 1 ]; then
    echo "Building STARS-H-CORE..."
    cd $BASEDIR/stars-h-core
    mkdir -p build && cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX
    make -j 4 -l 5|| make VERBOSE=1
    make install
fi

# HiCMA
echo "Building HiCMA..."
cd $BASEDIR/
rm -rf build
mkdir -p build && cd build
cmake .. -DCMAKE_C_FLAGS="-DMKL" -DHICMA_USE_MPI=$MPIVALUE -DBUILD_SHARED_LIBS=ON -DHICMA_ENABLE_TESTING=ON -DHICMA_ENABLE_TIMING=ON -DCMAKE_INSTALL_PREFIX=$PREFIX -DBLAS_LIBRARIES="-L${MKLROOT}/lib -Wl,-rpath,${MKLROOT}/lib -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"  $python
make -j 4 -l 5|| make VERBOSE=1
make install


cd $currdir
set +x
return


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

