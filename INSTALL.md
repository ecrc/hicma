# Installation on ECRC servers

1.  Go to home directory

2.  Get HiCMA from git repository

        git clone git@github.com:ecrc/hicma

3.  Run the following command: (This command will also setup PKG_CONFIG_PATH for STARS-H)

        . hicma/scripts/build.sh

4.  Run the following command for seeing which experimental cases will be executed:
    
        cd ./hicma
        ./exp/distruns.sh exp/cases/statistics.sh dry

5.  Remove the keyword "dry" to really run the cases.

# Installation

Installation requires `CMake` of version 3.2.3 at least. To build HiCMA,
follow these instructions:

1.  Get HiCMA from git repository

        git clone git@github.com:ecrc/hicma


2.  Go into hicma folder

        cd hicma

3.  Get submodules using git as follows. The submodules Chameleon, HCORE and STARS-H should be compiled and the `PKG_CONFIG_PATH` should be set.

        git submodule update --init --recursive

4.  Create build directory and go there

        mkdir build && cd build

5.  Use CMake to get all the dependencies

        cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install/  -DHICMA_USE_MPI=ON

6.  Build HiCMA

        make -j

7.  Build local documentation (optional)

        make docs

8.  Install HiCMA

        make install

9. Add line

        export PKG_CONFIG_PATH=/path/to/install:$PKG_CONFIG_PATH

    to your .bashrc file to use HiCMA as a library.

Now you can use `pkg-config` executable to collect compiler and linker flags for HiCMA or run the binaries under `/path/to/install/timing/`.

There are `scripts/build.sh` and `scripts/build-nompi.sh` scripts in the repository to build the whole software stack.

# CTests

You can run CTests in the build folder via the following command:

        ctest

A specific test can be run in verbose mode as follows:

        ctest -R time_zgemm_tile-mpi-edsin -V

# CMake Configuration

The following definitions affect how HiCMA library works:

**HICMA_ALWAYS_FIX_RANK**: Disables rank descriptors. The rank of the input matrix will be uniform 
across all tiles. This option can be used to measure the impact of additionally 
communicating rank descriptor during computations. If rank descriptor is not communicated, 
the total number of messages becomes less.

**HICMA_USE_MPI=[ON,OFF]**: Enables/Disables MPI in HiCMA.

