#!/bin/bash
repo=hicma-dev
hwloc_install_dir=$HOME/hicma-dev/hwloc-install-allinea
starpu_install_dir=$HOME/hicma-dev/starpu-1.2-install-allinea


currentdir=`pwd`
run_clone=0
set_pkgconfig_runtime_libs=0
run_update_submodules=0
run_module_setup=0
compile_hwloc=0
compile_starpu=0
compile_cham=0
compile_starsh=1
compile_starshcore=0
compile_hicma=1
pause_info(){
    echo "Please press enter key to proceed"; read
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
    . scripts/modules-isambard-allinea.sh
    pause_info
fi
if [ $compile_hwloc -eq 1 ];then
    if [  ! -f "hwloc-1.11.13.tar.gz" ]; then
	wget https://download.open-mpi.org/release/hwloc/v1.11/hwloc-1.11.13.tar.gz
        tar -zxvf hwloc-1.11.13.tar.gz
    fi
    if [ -d "hwloc-1.11.13" ]; then
    	rm -rf hwloc-1.11.13
    fi
    tar -zxvf hwloc-1.11.13.tar.gz
    cd hwloc-1.11.13
    [[ -d $hwloc_install_dir ]] || mkdir -p $hwloc_install_dir
    CC=cc CXX=CC ./configure --prefix=$hwloc_install_dir --disable-libxml2 --disable-pci --enable-shared=yes --enable-static=yes
    make -j
    make -j install
    if [ -d "$hwloc_install_dir/lib/pkgconfig" ]; then
        export PKG_CONFIG_PATH=$hwloc_install_dir/lib/pkgconfig:$PKG_CONFIG_PATH
    fi
    pause_info
fi
if [ $compile_starpu -eq 1 ];then
    if [ ! -f "starpu-1.2.6.tar.gz" ]; then
        wget http://starpu.gforge.inria.fr/files/starpu-1.2.6/starpu-1.2.6.tar.gz
    fi
    if [ -d "starpu-1.2.6" ]; then
        rm -rf starpu-1.2.6
    fi
    tar -zxvf starpu-1.2.6.tar.gz
    cd starpu-1.2.6
    [[ -d $starpu_install_dir ]] || mkdir -p $starpu_install_dir
    CC=cc CXX=CC ./configure --prefix=$starpu_install_dir --disable-cuda --disable-opencl  --disable-build-doc --disable-export-dynamic --with-mpicc=`which cc`
    #--disable-mpi-check 
    make -j
    make -j install
    if [ -d "$starpu_install_dir/lib/pkgconfig" ]; then
        export PKG_CONFIG_PATH=$starpu_install_dir/lib/pkgconfig:$PKG_CONFIG_PATH
    fi
    pause_info
fi
cd $currentdir/$repo
MKLROOT=$ARMPLROOT
MKLLIBS=$ARMPLLIBS
if [ $compile_cham -eq 1 ];then
    cd chameleon
    if [ -d build ]; then
        rm -rf build
    fi
    mkdir build 
    cd build
    cmake .. -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$currentdir/$repo/chameleon/build/install -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=OFF -DCHAMELEON_ENABLE_EXAMPLE=ON -DCHAMELEON_ENABLE_TESTING=OFF -DCHAMELEON_ENABLE_TIMING=ON -DHICMA_USE_MPI=ON -DCHAMELEON_USE_CUDA=OFF -DCHAMELEON_SCHED_QUARK=OFF -DCHAMELEON_SCHED_STARPU=ON -DSTARPU_DIR=$starpu_install_dir \
    -DBLAS_LIBRARIES="-Wl,--no-as-needed;-L${MKLROOT}/lib;${MKLLIBS};-lpthread;-lm;-ldl" -DBLAS_COMPILER_FLAGS="-I${MKLROOT}/include" -DLAPACK_LIBRARIES="-Wl,--no-as-needed;-L${MKLROOT}/lib;${MKLLIBS} -lpthread;-lm;-ldl" -DCBLAS_DIR="${MKLROOT}" -DLAPACKE_DIR="${MKLROOT}" -DTMG_DIR="${MKLROOT}" \
    -DMPI_C_COMPILER=`which cc`
    make -j install
    if [ -d $currentdir/$repo/chameleon/build/install ]; then
        export PKG_CONFIG_PATH=$currentdir/$repo/chameleon/build/install/lib/pkgconfig/:$PKG_CONFIG_PATH
    fi
    pause_info
fi
cd $currentdir/$repo
if [ $compile_starsh -eq 1 ];then
    cd stars-h
    rm -rf build
    mkdir build 
    cd build
    CFLAGS="-fPIC -fopenmp -mcpu=native -I${MKLROOT}/include" cmake .. -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$currentdir/$repo/stars-h/build/install -DMPI=OFF -DOPENMP=OFF -DSTARPU=OFF -DGSL=OFF -DBLAS_LIBRARIES="-L${MKLROOT}/lib ${MKLLIBS}"  -DBLAS_COMPILER_FLAGS="-fopenmp -mcpu=native -I${MKLROOT}/include" -DCBLAS_DIR="-L${MKLROOT}/lib ${MKLLIBS}"
    make -j install
    if [ -d $currentdir/$repo/stars-h/build/install ]; then
        export PKG_CONFIG_PATH=$currentdir/$repo/stars-h/build/install/lib/pkgconfig/:$PKG_CONFIG_PATH
    fi
    pause_info
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
    pause_info
fi
cd $currentdir/$repo
if [ $compile_hicma -eq 1 ];then
    mkdir exp/out
    mkdir exp/err
    rm -rf build
    mkdir build
    cd build
    CFLAGS="-DARMPL -I${MKLROOT}/include" cmake .. -DCMAKE_CXX_COMPILER=CC -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn -DCMAKE_INSTALL_PREFIX=$currentdir/$repo/build/install -DHICMA_USE_MPI=1 -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=OFF -DCHAMELEON_ENABLE_EXAMPLE=ON -DCHAMELEON_ENABLE_TESTING=ON -DCHAMELEON_ENABLE_TIMING=ON -DHICMA_USE_MPI=ON -DCHAMELEON_USE_CUDA=OFF -DCHAMELEON_SCHED_QUARK=OFF \
    -DCHAMELEON_SCHED_STARPU=ON -DSTARPU_DIR=$starpu_install_dir -DBLAS_LIBRARIES="-Wl,--no-as-needed;-L${MKLROOT}/lib;${MKLLIBS};-lpthread;-lm;-ldl" -DBLAS_COMPILER_FLAGS="-m64;-I${MKLROOT}/include" -DLAPACK_LIBRARIES="-Wl,--no-as-needed;-L${MKLROOT}/lib;${MKLLIBS};-lpthread;-lm;-ldl" -DCBLAS_DIR="${MKLROOT}" -DLAPACKE_DIR="${MKLROOT}" -DTMG_DIR="${MKLROOT}" \
    -DMPI_C_COMPILER=`which cc`
    make -j
fi

echo "Everything has finished, returning to initial folder \"$currentdir\""
cd $currentdir
