#!/bin/bash
repo=hicma-dev
hwloc_install_dir=/ibex/scratch/omairyrm/hicma-cs/hicma-dev/hwloc-install
starpu_install_dir=/ibex/scratch/omairyrm/hicma-cs/hicma-dev/starpu-1.2-install


currentdir=`pwd`
run_clone=0
set_pkgconfig_runtime_libs=1
run_update_submodules=0
run_module_setup=0
compile_hwloc=1
compile_starpu=1
compile_cham=1
compile_starsh=1
compile_starshcore=0
compile_hicma=1
compile_hcore=1
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
cd  $currentdir/$repo
if [ $run_update_submodules -eq 1 ];then 
    git submodule update --init --recursive
fi
if [ $run_module_setup -eq 1 ]; then
    . scripts/modules-isambard-allinea.sh
    pause_info
fi
cd  $currentdir/$repo
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
   ./configure --prefix=$hwloc_install_dir --disable-libxml2 --disable-pci --enable-shared=yes --enable-static=yes
    make -j
    make -j install
    if [ -d "$hwloc_install_dir/lib/pkgconfig" ]; then
        export PKG_CONFIG_PATH=$hwloc_install_dir/lib/pkgconfig:$PKG_CONFIG_PATH
    fi
    pause_info
fi
cd  $currentdir/$repo
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
    ./configure --prefix=$starpu_install_dir --disable-cuda --disable-opencl  --disable-build-doc --disable-export-dynamic  --without-mpicc --enable-maxcpus=40 
    #--disable-mpi-check 
    make -j
    make -j install
    if [ -d "$starpu_install_dir/lib/pkgconfig" ]; then
        export PKG_CONFIG_PATH=$starpu_install_dir/lib/pkgconfig:$PKG_CONFIG_PATH
    fi
    pause_info
fi
cd $currentdir/$repo
if [ $compile_cham -eq 1 ];then
    cd chameleon
    if [ -d build ]; then
        rm -rf build
    fi
    mkdir build 
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=$currentdir/$repo/chameleon/build/install -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON -DBUILD_SHARED_LIBS=OFF -DCHAMELEON_ENABLE_EXAMPLE=ON -DCHAMELEON_ENABLE_TESTING=OFF -DCHAMELEON_ENABLE_TIMING=ON -DHICMA_USE_MPI=OFF -DCHAMELEON_USE_CUDA=OFF -DCHAMELEON_SCHED_QUARK=OFF -DCHAMELEON_SCHED_STARPU=ON
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
    cmake .. -DCMAKE_INSTALL_PREFIX=$currentdir/$repo/stars-h/build/install -DMPI=OFF -DOPENMP=OFF -DSTARPU=OFF -DGSL=OFF
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
    cmake .. -DCMAKE_INSTALL_PREFIX=$currentdir/$repo/stars-h-core/build/install -DSTARSHCORE_WITH_GSL=OFF
    make -j install
    if [ -d $currentdir/$repo/stars-h-core/build/install ]; then
        export PKG_CONFIG_PATH=$currentdir/$repo/stars-h-core/build/install/lib/pkgconfig/:$PKG_CONFIG_PATH
    fi
    pause_info
fi
cd $currentdir/$repo
if [ $compile_hcore -eq 1 ];then
    cd hcore
    rm -rf build
    mkdir build
    cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=$PWD/install/
    make -j install
    if [ -d $currentdir/$repo/hcore/build/install ]; then
        export PKG_CONFIG_PATH=$currentdir/$repo/hcore/build/install/lib/pkgconfig/:$PKG_CONFIG_PATH
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
    cmake .. -DCMAKE_INSTALL_PREFIX=$currentdir/$repo/build/install  -DCMAKE_COLOR_MAKEFILE:BOOL=ON -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON
    make -j
fi

echo "Everything has finished, returning to initial folder \"$currentdir\""
cd $currentdir
