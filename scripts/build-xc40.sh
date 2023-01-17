#!/bin/bash
repo=hicma-dev
#starpu_install_dir=/project/k1205/akbudak/codes/starpu13-install-cdt
starpu_install_dir=/project/k1205/omairyrm/starpu-1.2.6/install
hwloc_install_dir=/project/k1205/omairyrm/hwloc/hwloc-1.11.10/install


currentdir=`pwd`
run_clone=0
set_pkgconfig_runtime_libs=1
run_update_submodules=0
run_module_setup=0
compile_cham=1
compile_starsh=1
compile_starshcore=0
compile_hicma=1
pause_info(){
    #echo "Please press enter key to proceed"; read
    echo "================="
}
pause_step(){
    echo "Please press enter key to proceed"; read
    #echo "================="
}

if [ $run_clone -eq 1 ]; then
if [ ! -d $repo ]; then
    git clone git@github.com:ecrc/$repo.git
    pause_step
else
    echo "Echo \"$repo\" exists so I am not cloning from github"
    pause_info
fi
fi
if [ $set_pkgconfig_runtime_libs -eq 1 ]; then
    export PKG_CONFIG_PATH=$hwloc_install_dir/lib/pkgconfig:$PKG_CONFIG_PATH
    export PKG_CONFIG_PATH=$starpu_install_dir/lib/pkgconfig:$PKG_CONFIG_PATH
fi
cd $repo
if [ $run_update_submodules -eq 1 ];then 
    git submodule update --init --recursive
fi
if [ $run_module_setup -eq 1 ]; then
    . scripts/modules-xc40.sh
fi
if [ $compile_cham -eq 1 ];then
    cd chameleon
    mkdir build 
    cd build
    cmake .. -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$currentdir/$repo/chameleon/build/install -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=OFF -DCHAMELEON_ENABLE_EXAMPLE=ON -DCHAMELEON_ENABLE_TESTING=ON -DCHAMELEON_ENABLE_TIMING=ON -DHICMA_USE_MPI=ON -DCHAMELEON_USE_CUDA=OFF -DCHAMELEON_SCHED_QUARK=OFF -DCHAMELEON_SCHED_STARPU=ON -DSTARPU_DIR=$starpu_install_dir \
    -DBLAS_LIBRARIES="-Wl,--no-as-needed;-L${MKLROOT}/lib;-lmkl_intel_lp64;-lmkl_core;-lmkl_sequential;-lpthread;-lm;-ldl" -DBLAS_COMPILER_FLAGS="-m64;-I${MKLROOT}/include" -DLAPACK_LIBRARIES="-Wl,--no-as-needed;-L${MKLROOT}/lib;-lmkl_intel_lp64;-lmkl_core;-lmkl_sequential;-lpthread;-lm;-ldl" -DCBLAS_DIR="${MKLROOT}" -DLAPACKE_DIR="${MKLROOT}" -DTMG_DIR="${MKLROOT}" \
    -DMPI_C_COMPILER=`which cc`
    make -j install
    if [ -d $currentdir/$repo/chameleon/build/install ]; then
        export PKG_CONFIG_PATH=$currentdir/$repo/chameleon/build/install/lib/pkgconfig/:$PKG_CONFIG_PATH
    fi
fi
cd $currentdir/$repo
if [ $compile_starsh -eq 1 ];then
    cd stars-h
    rm -rf build
    mkdir build 
    cd build
    CFLAGS=-fPIC cmake .. -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$currentdir/$repo/stars-h/build/install -DMPI=OFF -DOPENMP=OFF -DSTARPU=OFF -DGSL=OFF
    make -j install
    if [ -d $currentdir/$repo/stars-h/build/install ]; then
        export PKG_CONFIG_PATH=$currentdir/$repo/stars-h/build/install/lib/pkgconfig/:$PKG_CONFIG_PATH
    fi
fi
cd $currentdir/$repo
if [ $compile_starshcore -eq 1 ];then
    cd stars-h-core
    rm -rf build
    mkdir build 
    cd build
    CFLAGS=-fPIC cmake .. -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$currentdir/$repo/stars-h-core/build/install -DSTARSHCORE_WITH_GSL=OFF
    make -j install
    if [ -d $currentdir/$repo/stars-h-core/build/install ]; then
        export PKG_CONFIG_PATH=$currentdir/$repo/stars-h-core/build/install/lib/pkgconfig/:$PKG_CONFIG_PATH
    fi
fi

# HCORE
cd $currentdir/$repo
cd hcore
rm -rf build
mkdir -p build/installdir
cd build
   CFLAGS=-fPIC cmake ..  -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn  -DCMAKE_INSTALL_PREFIX=$PWD/installdir 
make clean
make -j
make install
export PKG_CONFIG_PATH=$currentdir/$repo/hcore/build/installdir/lib/pkgconfig:$PKG_CONFIG_PATH


cd $currentdir/$repo
if [ $compile_hicma -eq 1 ];then
    mkdir exp/out
    mkdir exp/err
    rm -rf build
    mkdir build
    cd build
    cmake .. -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$currentdir/$repo/build/install -DHICMA_USE_MPI=1 -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=OFF -DCHAMELEON_ENABLE_EXAMPLE=ON -DCHAMELEON_ENABLE_TESTING=ON -DCHAMELEON_ENABLE_TIMING=ON -DHICMA_USE_MPI=ON -DCHAMELEON_USE_CUDA=OFF -DCHAMELEON_SCHED_QUARK=OFF \
    -DCHAMELEON_SCHED_STARPU=ON -DSTARPU_DIR=$starpu_install_dir -DBLAS_LIBRARIES="-Wl,--no-as-needed;-L${MKLROOT}/lib;-lmkl_intel_lp64;-lmkl_core;-lmkl_sequential;-lpthread;-lm;-ldl" -DBLAS_COMPILER_FLAGS="-m64;-I${MKLROOT}/include" -DLAPACK_LIBRARIES="-Wl,--no-as-needed;-L${MKLROOT}/lib;-lmkl_intel_lp64;-lmkl_core;-lmkl_sequential;-lpthread;-lm;-ldl" -DCBLAS_DIR="${MKLROOT}" -DLAPACKE_DIR="${MKLROOT}" -DTMG_DIR="${MKLROOT}" \
    -DMPI_C_COMPILER=`which cc`
    make -j
fi

echo "Everything has finished, returning to initial folder \"$currentdir\""
cd $currentdir
