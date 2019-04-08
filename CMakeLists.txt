###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2016 Inria. All rights reserved.
# @copyright (c) 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
# @copyright (c) 2017      King Abdullah University of Science and Technology (KAUST)
#
###
#
#  @file CMakeLists.txt
#
#  @project HiCMA
#  HiCMA is a software package provided by:
#     King Abdullah University of Science and Technology
#
#  @version 0.1.1
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @author Eduardo Gonzalez
#  @author Kadir Akbudak
#  @date 2018-11-08
#
## disables all rank descriptors and runs in fixed rank mode
#add_definitions(-DHICMA_ALWAYS_FIX_RANK)

## disables all computations so affect of communication can be measured
#add_definitions(-DHICMA_DISABLE_ALL_COMPUTATIONS)

## disables low rank operations except the dense potrf
#add_definitions(-DHICMA_DISABLE_HCORE_COMPUTATIONS)

#LAPACKE_dge_trans() is used in HiCMA. Appropriate header file must be included
#add_definitions(-DMKL)
#add_definitions(-DLAPACKE_UTILS)

cmake_minimum_required(VERSION 3.2.3)

# directly make an error if in-source build
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
   message(FATAL_ERROR "In-source builds are not allowed.\n"
   "Please create a build directory first and execute cmake configuration from "
   "this directory. Example: mkdir build && cd build && cmake ..")
endif()

project(HICMA C )

# set project version number
set(HICMA_VERSION_MAJOR 0)
set(HICMA_VERSION_MINOR 1)
set(HICMA_VERSION_PATCH 1)

set(MORSE_CMAKE_DIR "" CACHE PATH "Directory of MORSE CMake modules, can be external to the project")

## CMAKE MODULES :: ECRC
## REQUIRED FOR TESTS TO LINK LIBRARIES
if(NOT EXISTS "${PROJECT_SOURCE_DIR}/cmake_modules/ecrc/modules" )
    find_package(Git REQUIRED)
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule init WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_init OUTPUT_QUIET ERROR_QUIET)
    execute_process(COMMAND ${GIT_EXECUTABLE} submodule update WORKING_DIRECTORY ${PROJECT_SOURCE_DIR} RESULT_VARIABLE _res_update OUTPUT_QUIET ERROR_QUIET)
    if( ${_res_init} GREATER 0 OR ${_res_update} GREATER 0 )
        message(FATAL_ERROR "ECRC CMake modules were not found.\n"
                            "We tried: 'git submodule init && git submodule update' and resulted in error" )
    endif()
endif()
## ECRC INITIALIZATION
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules)
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake_modules/ecrc/modules")
set(ECRC_CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake_modules/ecrc/modules )
include(EcrcInit)
include(GenPkgConfig)


# # Add extra cmake module path and initialize morse cmake modules
# # --------------------------------------------------------------
# list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules)
# if(MORSE_CMAKE_DIR)
#     list(APPEND CMAKE_MODULE_PATH "${MORSE_CMAKE_DIR}/cmake_modules/morse_cmake/modules")
#     set(MORSE_CMAKE_MODULE_PATH ${MORSE_CMAKE_DIR}/cmake_modules/morse_cmake/modules )
# elseif(EXISTS "${CMAKE_SOURCE_DIR}/cmake_modules/morse_cmake")
#     list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake_modules/morse_cmake/modules")
#     set(MORSE_CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake_modules/morse_cmake/modules )
# else()
#   message(FATAL_ERROR "HiCMA CMake system relies on morse_cmake modules developed here: "
#     "https://gitlab.inria.fr/solverstack/morse_cmake. Please set MORSE_CMAKE_DIR to this source "
#     "directory.")
# endif()
# include(MorseInit)
# include(GenPkgConfig)

#############################################
#                                           #
#        Compilation of HICMA           #
#                                           #
#############################################


###############################################################################
# Parameters/Options #
######################

# Set the RPATH config
# --------------------

# use, i.e. don't skip the full RPATH for the build tree
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# when building, don't use the install RPATH already
# (but later on when installing)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

# the RPATH to be used when installing
list(APPEND CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

# Misc options
# ------------
option(BUILD_SHARED_LIBS "Build shared libraries" OFF)

# Define precision supported by HICMA
# -----------------------------------------
set( RP_HICMA_DICTIONNARY ${PROJECT_SOURCE_DIR}/cmake_modules/morse_cmake/modules/precision_generator/subs.py )
set( RP_HICMA_PRECISIONS  "s;d;c;z" )
include(RulesPrecisions)

# Options to select the runtime
# -----------------------------

# Create a list of possible runtime
#set(HICMA_SCHED_list "PARSEC" "STARPU" "QUARK")
set(HICMA_SCHED_list "STARPU")

# Initially PaRSEC runtime is enabled
option(HICMA_SCHED_PARSEC
    "Enable PaRSEC scheduler as the default runtime
    (Conflict with other HICMA_SCHED_* options)" OFF)
option(HICMA_SCHED_STARPU
    "Enable StarPU scheduler as the default runtime
    (Conflict with other HICMA_SCHED_* options)" ON)
option(HICMA_SCHED_QUARK
    "Enable Quark scheduler as the default runtime
    (Conflict with other HICMA_SCHED_* options)" OFF)

# For now, we are able to compile only one runtime at a time, so we disable combinations of runtimes
if (HICMA_SCHED_STARPU)
    set(HICMA_SCHED_QUARK OFF)
    set(HICMA_SCHED_PARSEC OFF)
elseif (HICMA_SCHED_QUARK)
    set(HICMA_SCHED_STARPU OFF)
    set(HICMA_SCHED_PARSEC OFF)
elseif (HICMA_SCHED_PARSEC)
    set(HICMA_SCHED_QUARK OFF)
    set(HICMA_SCHED_STARPU OFF)
endif()
# Set default to StarPU if nothing specific is required by the user
if ( NOT HICMA_SCHED_STARPU AND NOT HICMA_SCHED_PARSEC AND NOT HICMA_SCHED_QUARK )
    set(HICMA_SCHED_STARPU ON)
endif()
if (HICMA_SCHED_STARPU)
    message("-- ${BoldGreen}HICMA_SCHED_STARPU is set to ON: HICMA uses StarPU runtime\n"
            "   To use HICMA with Quark  runtime: set HICMA_SCHED_QUARK  to ON\n"
            "   To use HICMA with PaRSEC runtime: set HICMA_SCHED_PARSEC to ON\n"
            "   (HICMA_SCHED_STARPU will be disabled)${ColourReset}")
elseif(HICMA_SCHED_QUARK)
    message("-- ${BoldGreen}HICMA_SCHED_QUARK is set to ON: HICMA uses Quark runtime\n"
            "   To use HICMA with StarPU runtime: set HICMA_SCHED_STARPU to ON\n"
            "   To use HICMA with PaRSEC runtime: set HICMA_SCHED_PARSEC to ON\n"
            "   (HICMA_SCHED_QUARK will be disabled)${ColourReset}")
elseif(HICMA_SCHED_PARSEC)
    message("-- ${BoldGreen}HICMA_SCHED_PARSEC is set to ON: HICMA uses PaRSEC runtime\n"
            "   To use HICMA with StarPU runtime: set HICMA_SCHED_STARPU to ON\n"
            "   To use HICMA with Quark  runtime: set HICMA_SCHED_QUARK  to ON\n"
            "   (HICMA_SCHED_PARSEC will be disabled)${ColourReset}")
endif()

# Check that one, and only one, SCHED option is set to ON
# count number of runtime sets to ON
math(EXPR number_of_active_runtime 0)
foreach (runtime ${HICMA_SCHED_list})
    if (HICMA_SCHED_${runtime})
        math(EXPR number_of_active_runtime "${number_of_active_runtime}+1")
    endif()
endforeach()
if (NOT number_of_active_runtime STREQUAL 1)
    message(FATAL_ERROR
            "Number of active runtime is ${number_of_active_runtime}, "
            "the user should activate one (and only one) runtime. ")
endif()

# Use intermediate variable since cmake_dependent_option doesn't have OR conditions
set(HICMA_ENABLE_MPI OFF CACHE INTERNAL "Tells if MPI might be supported by the runtime")
if ( HICMA_SCHED_PARSEC OR HICMA_SCHED_STARPU )
     set(HICMA_ENABLE_MPI ON FORCE)
endif()

# # Use intermediate variable since cmake_dependent_option doesn't have OR conditions
# set(HICMA_ENABLE_CUDA OFF CACHE INTERNAL "Tells if CUDA might be supported by the runtime")
# if ( HICMA_SCHED_PARSEC OR HICMA_SCHED_STARPU )
#      set(HICMA_ENABLE_CUDA ON FORCE)
# endif()

# Additional options
# ------------------

# Enable the distributed interface (allowed only when StarPU or PaRSEC is enabled)
# TODO: Default should be changed to ON/OFF when it will be ok
cmake_dependent_option(HICMA_USE_MPI
                    "Enable distributed memory through MPI" OFF
                    "HICMA_ENABLE_MPI" OFF)
if (HICMA_ENABLE_MPI AND NOT HICMA_USE_MPI)
    message("-- ${BoldGreen}HICMA_USE_MPI is set to OFF, turn it ON to use MPI (unsupported by Quark)${ColourReset}")
endif()

# # Enable Cuda kernels if StarPU (only if StarPU or PaRSEC runtime is enabled)
# # TODO: Default should be changed to ON/OFF when it will be ok
# cmake_dependent_option(HICMA_USE_CUDA
#                     "Enable CUDA kernels" OFF
#                     "HICMA_ENABLE_CUDA" OFF)
# if (HICMA_ENABLE_CUDA AND NOT HICMA_USE_CUDA)
#     message("-- ${BoldGreen}HICMA_USE_CUDA is set to OFF, turn it ON to use CUDA (unsupported by Quark)${ColourReset}")
# endif()

# # Enable MAGMA advanced kernels if CUDA is enabled
# cmake_dependent_option(HICMA_USE_MAGMA
#                     "Enable MAGMA Cuda kernels" OFF
#                     "HICMA_USE_CUDA" OFF)
# if (HICMA_USE_CUDA AND NOT HICMA_USE_MAGMA)
#     message("-- ${BoldGreen}HICMA_USE_MAGMA is set to OFF, turn it ON to use MAGMA (only with StarPU)${ColourReset}")
# endif()

# Enable FXT if StarPU
option(HICMA_ENABLE_TRACING "Enable tracing support" OFF)
if (NOT HICMA_ENABLE_TRACING)
    message("-- ${BoldGreen}HICMA_ENABLE_TRACING is set to OFF, turn it ON to use FxT (with StarPU)${ColourReset}")
endif()
#option(HICMA_USE_EZTRACE "Enable EZTRACE to build modules" OFF)

# # Options to enable/disable testings and timings
# # ----------------------------------------------
# option(HICMA_ENABLE_DOCS     "Enable documentation build"  OFF)
# if (HICMA_ENABLE_DOCS)
#     message("-- ${BoldGreen}HICMA_ENABLE_DOCS is set to ON, turn it OFF to avoid building docs${ColourReset}")
# endif()
# option(HICMA_ENABLE_EXAMPLE  "Enable examples build"       ON)
# if (HICMA_ENABLE_EXAMPLE)
#     message("-- ${BoldGreen}HICMA_ENABLE_EXAMPLE is set to ON, turn it OFF to avoid building examples${ColourReset}")
# endif()
# option(HICMA_ENABLE_TESTING  "Enable testings build"       ON)
# if (HICMA_ENABLE_TESTING)
#     message("-- ${BoldGreen}HICMA_ENABLE_TESTING is set to ON, turn it OFF to avoid building testing${ColourReset}")
# endif()
option(HICMA_ENABLE_DOCS "Build documentation in docs directory" ON)
# If docs option is ON
if(HICMA_ENABLE_DOCS)
    find_package(Doxygen)
    if(DOXYGEN_FOUND)
        add_subdirectory("docs")
    else()
        message(STATUS "Doxygen NOT found, skipping it")
    endif()
endif()
option(HICMA_ENABLE_TIMING   "Enable timings build"        ON)
if (HICMA_ENABLE_TIMING)
    message("-- ${BoldGreen}HICMA_ENABLE_TIMING is set to ON, turn it OFF to avoid building timing${ColourReset}")
endif()
option(HICMA_ENABLE_TESTING   "Enable testings build"        ON)
if (HICMA_ENABLE_TESTING)
    message("-- ${BoldGreen}HICMA_ENABLE_TESTING is set to ON, turn it OFF to avoid building testing${ColourReset}")
endif()

# Option to activate or not simulation mode (use Simgrid through StarPU)
# ----------------------------------------------------------------------
if(HICMA_SCHED_STARPU)
    option(HICMA_SIMULATION "Enable simulation mode using Simgrid through StarPU" OFF)
    if (NOT HICMA_SIMULATION)
        message("-- ${BoldGreen}HICMA_SIMULATION is set to OFF, turn it ON to use"
            " SIMULATION mode (only with StarPU compiled with SimGrid)${ColourReset}")
    endif()
    option(HICMA_ENABLE_PRUNING_STATS "Enable pruning statistics" OFF)
    if (NOT HICMA_ENABLE_PRUNING_STATS)
        message("-- ${BoldGreen}HICMA_ENABLE_PRUNING_STATS is set to OFF, turn it ON to build pruning statistics${ColourReset}")
    endif()
    if(HICMA_ENABLE_PRUNING_STATS)
        add_definitions(-DHICMA_ENABLE_PRUNING_STATS)
    endif(HICMA_ENABLE_PRUNING_STATS)
endif()

# Initially we need to generate files for different precisions
# TODO: use this option to avoid generate multiple precisions each time we launch cmake
#option(HICMA_GEN_PREC "Generate source files precisions" ON)
#------------------------------------------------------------------------------


###############################################################################
# Look for dependencies #
#########################
set(HICMA_DEP "")

# Check for Thread library
# ------------------------
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package(Threads REQUIRED)
if( THREADS_FOUND )
    list(APPEND EXTRA_LIBRARIES ${CMAKE_THREAD_LIBS_INIT})
endif ()

# Add math library to the list of extra
# it normally exists on all common systems provided with a C compiler
set(M_LIBRARIES "")
if(UNIX OR WIN32)
    find_library(
        M_m_LIBRARY
        NAMES m
        )
    mark_as_advanced(M_m_LIBRARY)
    if (M_m_LIBRARY)
        list(APPEND M_LIBRARIES "${M_m_LIBRARY}")
        list(APPEND EXTRA_LIBRARIES "${M_m_LIBRARY}")
    else()
        message(FATAL_ERROR "Could NOT find libm on your system."
            " Are you sure to a have a C compiler installed?")
    endif()
endif()

# Try to find librt (libposix4 - POSIX.1b Realtime Extensions library)
# on Unix systems except Apple ones because it does not exist on it
set(RT_LIBRARIES "")
if(UNIX AND NOT APPLE)
    find_library(
        RT_rt_LIBRARY
        NAMES rt
        )
    mark_as_advanced(RT_rt_LIBRARY)
    if (RT_rt_LIBRARY)
        list(APPEND RT_LIBRARIES "${RT_rt_LIBRARY}")
        list(APPEND EXTRA_LIBRARIES "${RT_rt_LIBRARY}")
    else()
        message(FATAL_ERROR "Could NOT find librt on your system")
    endif()
endif()

# If simulation we don't enter in kernel functions so that we don't need to
# link with concerned libraries
if(NOT HICMA_SIMULATION)

    # HICMA depends on CBLAS
    #---------------------------
    find_package(CBLAS COMPONENTS BLASEXT)
    if(BLAS_FOUND)
        if (BLAS_LIBRARY_DIRS)
            # the RPATH to be used when installing
            list(APPEND CMAKE_INSTALL_RPATH "${BLAS_LIBRARY_DIRS}")
        endif()
        if(BLAS_LINKER_FLAGS)
            list(APPEND CMAKE_EXE_LINKER_FLAGS "${BLAS_LINKER_FLAGS}")
        endif()
    else()
        message(FATAL_ERROR "BLAS library has not been found")
    endif()
    if(CBLAS_FOUND)
        include_directories(${CBLAS_INCLUDE_DIRS})
        if(CBLAS_LIBRARY_DIRS)
            # the RPATH to be used when installing
            list(APPEND CMAKE_INSTALL_RPATH "${CBLAS_LIBRARY_DIRS}")
        endif()
        if (CBLAS_LIBRARIES)
            if (CBLAS_LIBRARIES_DEP)
                list(INSERT HICMA_DEP 0 ${CBLAS_LIBRARIES_DEP})
            else()
                list(INSERT HICMA_DEP 0 ${CBLAS_LIBRARIES})
            endif()
        endif()
    else()
        if(ECRC_VERBOSE_FIND_PACKAGE)
            if(CBLAS_STANDALONE OR NOT CBLAS_WORKS)
                if (NOT CBLAS_cblas.h_DIRS)
                    Print_Find_Header_Status(cblas cblas.h)
                endif ()
                if (NOT CBLAS_cblas_LIBRARY)
                    Print_Find_Library_Status(cblas libcblas)
                endif ()
            endif()
        else()
            message(WARNING "CBLAS library has not been found and ECRC_VERBOSE_FIND_PACKAGE is set to OFF."
            " Try to activate ECRC_VERBOSE_FIND_PACKAGE option (-DECRC_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
        endif()
        message(FATAL_ERROR "A CBLAS library is required but has not been found")
    endif()

    list(REMOVE_DUPLICATES CMAKE_EXE_LINKER_FLAGS)
    string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

    # HICMA depends on LAPACKE
    #-----------------------------
    # standalone version of lapacke seems useless for now
    # let the comment in case we meet some problems of non existing lapacke
    # functions in lapack library such as mkl, acml, ...
    #set(LAPACKE_STANDALONE TRUE)
    find_package(LAPACKE COMPONENTS LAPACKEXT)
    if(LAPACK_FOUND AND LAPACK_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${LAPACK_LIBRARY_DIRS}")
    else()
        message(FATAL_ERROR "A LAPACK library is required but has not been found")
    endif()
    if(LAPACKE_FOUND)
        include_directories(${LAPACKE_INCLUDE_DIRS})
        if(LAPACKE_LIBRARY_DIRS)
            # the RPATH to be used when installing
            list(APPEND CMAKE_INSTALL_RPATH "${LAPACKE_LIBRARY_DIRS}")
        endif()
        if(LAPACKE_LINKER_FLAGS)
            list(APPEND CMAKE_EXE_LINKER_FLAGS "${LAPACKE_LINKER_FLAGS}")
        endif()
        if (LAPACKE_LIBRARIES)
            if (LAPACKE_LIBRARIES_DEP)
                list(INSERT HICMA_DEP 0 ${LAPACKE_LIBRARIES_DEP})
            else()
                list(INSERT HICMA_DEP 0 ${LAPACKE_LIBRARIES})
            endif()
        endif()
    else()
        if(ECRC_VERBOSE_FIND_PACKAGE)
            if (LAPACKE_STANDALONE OR NOT LAPACKE_WORKS)
                if (NOT LAPACKE_lapacke.h_DIRS)
                    Print_Find_Header_Status(lapacke lapacke.h)
                endif ()
                if (NOT LAPACKE_lapacke_LIBRARY)
                    Print_Find_Library_Status(lapacke liblapacke)
                endif ()
            endif()
        else()
            message(WARNING "LAPACKE library has not been found and ECRC_VERBOSE_FIND_PACKAGE is set to OFF."
            " Try to activate ECRC_VERBOSE_FIND_PACKAGE option (-DECRC_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
        endif()
        message(FATAL_ERROR "A LAPACKE library is required but has not been found")
    endif()

#     # HICMA depends on TMG
#     #-------------------------
#     if (HICMA_ENABLE_TESTING OR HICMA_ENABLE_TIMING)
#         find_package(TMG)
#         if(TMG_FOUND AND TMG_LIBRARY_DIRS)
#             # the RPATH to be used when installing
#             list(APPEND CMAKE_INSTALL_RPATH "${TMG_LIBRARY_DIRS}")
#         else()
#             message(FATAL_ERROR "A TMG library is required but has not been found")
#         endif()
#         if(TMG_LINKER_FLAGS)
#             list(APPEND CMAKE_EXE_LINKER_FLAGS "${TMG_LINKER_FLAGS}")
#         endif()
#         if (TMG_LIBRARIES)
#             if (TMG_LIBRARIES_DEP)
#                 list(INSERT HICMA_DEP 0 ${TMG_LIBRARIES_DEP})
#             else()
#                 list(INSERT HICMA_DEP 0 ${TMG_LIBRARIES})
#             endif()
#         endif()
#     endif()

#     # HICMA depends on CUDA/CUBLAS
#     #----------------------------------
#     if (HICMA_USE_CUDA)
#
#         find_package(CUDA REQUIRED)
#
#         if (CUDA_FOUND)
#             message("-- ${Blue}Add definition HICMA_USE_CUDA"
#             " - Activate CUDA in HiCMA${ColourReset}")
#             set(EXTRA_LIBRARIES_CUDA)
#             set(HICMA_USE_CUDA 1)
#             if(CUDA_VERSION VERSION_LESS "3.0")
#                 set(CUDA_HOST_COMPILATION_CPP OFF)
#             endif(CUDA_VERSION VERSION_LESS "3.0")
#             set(CUDA_BUILD_EMULATION OFF)
#             if (CUDA_INCLUDE_DIRS)
#                 include_directories(${CUDA_INCLUDE_DIRS})
#             else()
#                 message(WARNING "HICMA_USE_CUDA requires"
#                 "\n   CUDA_INCLUDE_DIRS to be found. Be sure you have"
#                 "\n   cuda headers with your distribution of CUDA.")
#             endif()
#             if (CUDA_LIBRARIES)
#                 list(APPEND EXTRA_SOURCES src/gpu_data.c)
#                 set(saved_include "${CMAKE_REQUIRED_INCLUDES}")
#                 set(saved_libs "${CMAKE_REQUIRED_LIBRARIES}")
#                 set(CMAKE_REQUIRED_INCLUDES  "${CMAKE_REQUIRED_INCLUDES};${CUDA_INCLUDE_DIRS}")
#                 set(CMAKE_REQUIRED_LIBRARIES "${CMAKE_REQUIRED_LIBRARIES};${CUDA_LIBRARIES}")
#                 if(CUDA_VERSION VERSION_LESS "4.0")
#                     set(CUDA_HAVE_PEER_DEVICE_MEMORY_ACCESS 0)
#                 else()
#                     check_function_exists(cuDeviceCanAccessPeer CUDA_HAVE_PEER_DEVICE_MEMORY_ACCESS)
#                 endif()
#                 set(CMAKE_REQUIRED_INCLUDES  "${saved_include}")
#                 set(CMAKE_REQUIRED_LIBRARIES "${saved_libs}")
#                 # Add cublas if found
#                 if (CUDA_CUBLAS_LIBRARIES)
#                     message("-- ${Blue}Add definition HICMA_USE_CUBLAS"
#                     " - Use GPU kernels from cuBLAS${ColourReset}")
#                     set( HICMA_USE_CUBLAS 1)
#                     list(APPEND EXTRA_LIBRARIES_CUDA ${CUDA_CUBLAS_LIBRARIES})
#                     list(APPEND EXTRA_LIBRARIES_CUDA ${CUDA_LIBRARIES}) # need CUDA libs to link.
#                 else()
#                     message(FATAL_ERROR "HICMA_USE_CUDA requires"
#                     "\n   CUDA_CUBLAS_LIBRARIES to be found. Be sure you have"
#                     "\n   libcublas with your distribution of CUDA.")
#                 endif()
#             else()
#                 message(FATAL_ERROR "HICMA_USE_CUDA requires"
#                 "\n   CUDA_LIBRARIES to be found. Be sure you have"
#                 "\n   libcuda with your distribution of CUDA.")
#             endif()
#             #FindHeader(CUDA cublas_v2.h)
#             #if(CUDA_cublas_v2.h_DIRS)
#                 #set(HICMA_USE_CUBLAS_V2 1)
#             #endif()
#         endif (CUDA_FOUND)
#
#         # HICMA depends on MAGMA gpu kernels
#         # call our cmake module to test (in cmake_modules)
#         # change this call position if not appropriated
#         #----------------------------------------------
#         if ( CUDA_FOUND AND HICMA_USE_MAGMA )
#             set(HICMA_MAGMA_VERSION "1.4" CACHE STRING "oldest MAGMA version desired")
#             find_package(MAGMA ${HICMA_MAGMA_VERSION})
#             if ( MAGMA_FOUND )
#                 message("-- ${Blue}Add definition HICMA_USE_MAGMA"
#                 " - Use GPU kernels from MAGMA${ColourReset}")
#                 set(HICMA_USE_MAGMA 1)
#                 if(MAGMA_INCLUDE_DIRS)
#                     include_directories(${MAGMA_INCLUDE_DIRS})
#                 endif()
#                 if(MAGMA_LIBRARY_DIRS)
#                     # the RPATH to be used when installing
#                     list(APPEND CMAKE_INSTALL_RPATH "${MAGMA_LIBRARY_DIRS}")
#                 endif()
#                 if (MAGMA_LIBRARIES)
#                     list(INSERT EXTRA_LIBRARIES_CUDA 0 ${MAGMA_LIBRARIES})
#                 endif()
#             else( MAGMA_FOUND )
#                 if(ECRC_VERBOSE_FIND_PACKAGE)
#                     if (NOT MAGMA_magma.h_DIRS)
#                         Print_Find_Header_Status(magma magma.h)
#                     endif ()
#                     if (NOT MAGMA_magma_LIBRARY)
#                         Print_Find_Library_Status(magma libmagma)
#                     endif ()
#                 else()
#                     message(WARNING "MAGMA library has not been found and ECRC_VERBOSE_FIND_PACKAGE is set to OFF."
#                     " Try to activate ECRC_VERBOSE_FIND_PACKAGE option (-DECRC_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
#                 endif()
#                 message(FATAL_ERROR "MAGMA library is required but has not been found")
#             endif( MAGMA_FOUND )
#         endif()
#
#         list(APPEND EXTRA_LIBRARIES ${EXTRA_LIBRARIES_CUDA})
#
#     endif(HICMA_USE_CUDA)

    # HICMA depends on MPI
    #-------------------------
    if (HICMA_USE_MPI)

        # allows to use an external mpi compilation by setting compilers with
        # -DMPI_C_COMPILER=path/to/mpicc -DMPI_Fortran_COMPILER=path/to/mpif90
        # at cmake configure
        if(NOT MPI_C_COMPILER)
            set(MPI_C_COMPILER mpicc)
        endif()
        find_package(MPI REQUIRED)

        if (MPI_C_FOUND)
            message("-- ${Blue}Add definition HICMA_USE_MPI"
            " - Activate MPI in HiCMA${ColourReset}")
            set(HICMA_USE_MPI 1)
            list(APPEND EXTRA_LIBRARIES ${MPI_C_LIBRARIES} )
            include_directories( ${MPI_C_INCLUDE_PATH} )
            # tests for intel mpi
            #list(APPEND MPI_C_COMPILE_FLAGS "-mt_mpi")
            #list(APPEND MPI_COMPILE_FLAGS "-mt_mpi")
            if(MPI_C_LINK_FLAGS)
                list(APPEND CMAKE_EXE_LINKER_FLAGS "${MPI_C_LINK_FLAGS}")
            endif()
        endif (MPI_C_FOUND)

    endif (HICMA_USE_MPI)

else (NOT HICMA_SIMULATION)

    message("-- ${BoldBlue}Simulation mode of HICMA is activated (HICMA_SIMULATION=ON)."
    "\n   This mode allows you to simulate execution of algorithms with StarPU compiled with SimGrid.${ColourReset}")

    message("-- ${Blue}Add definition HICMA_SIMULATION"
    " - Activate simulation mode in HiCMA (to use with StarPU+SimGrid)${ColourReset}")
    set(HICMA_SIMULATION 1)

    add_subdirectory(simucore)

    # Guard against mis-used simulation mode
#     if(NOT DEFINED HICMA_USE_CUDA)
#         message(FATAL_ERROR "${BoldBlue}"
#         "In simulation mode HICMA_USE_CUDA and HICMA_USE_MAGMA should be set to"
#         "   ensure that the user is aware of the version to be used. If HiCMA's"
#         "   kernels are available for NVIDIA CUDA GPUs and if the according"
#         "   perfmodels are available in STARPU_HOME then use HICMA_USE_CUDA=ON"
#         "   else set HICMA_USE_CUDA=OFF. The same idea is applicable with MAGMA."
#         "${ColourReset}")
#     endif()
#     if(NOT DEFINED HICMA_USE_MAGMA)
#         message(WARNING "${BoldBlue}"
#         "In simulation mode HICMA_USE_CUDA and HICMA_USE_MAGMA should be set to"
#         "   ensure that the user is aware of the version to be used. If HiCMA's"
#         "   MAGMA kernels are available for NVIDIA CUDA GPUs and if the according"
#         "   perfmodels are available in STARPU_HOME then use HICMA_USE_MAGMA=ON"
#         "   else set HICMA_USE_MAGMA=OFF.${ColourReset}")
#     endif()

#     # Add CUDA definition if required
#     if (HICMA_USE_CUDA)
#         set(HICMA_USE_CUDA 1)
#         # Add MAGMA definition if required
#         if (HICMA_USE_MAGMA)
#             set(HICMA_USE_MAGMA 1)
#         endif()
#     endif()

    if (NOT HICMA_SCHED_STARPU)
        message(FATAL_ERROR "Simulation mode of HiCMA only works with"
        "\n   the StarPU runtime. Please use HICMA_SCHED_STARPU=ON.")
    endif()

    if (HICMA_ENABLE_EXAMPLE)
        set(HICMA_ENABLE_EXAMPLE OFF)
        message("-- ${BoldBlue}HICMA_ENABLE_EXAMPLE is set to ON but we turn it OFF.")
    endif ()
    if (HICMA_ENABLE_TESTING)
        set(HICMA_ENABLE_TESTING OFF)
        message("-- ${BoldBlue}HICMA_ENABLE_TESTING is set to ON but we turn it OFF."
            "\n   Because we are compiling the simulation mode (HICMA_SIMULATION=ON),"
            "\n   there is no sense in compiling testing drivers that are used to check"
            "\n   numerical correctness of algorithms and kernels.${ColourReset}")
    endif ()

    # Simulation mode: we depend on SimGrid
    find_package(SIMGRID REQUIRED)

    # HICMA depends on MPI
    #-------------------------
    if (HICMA_USE_MPI)

        # allows to use an external mpi compilation by setting compilers with
        # -DMPI_C_COMPILER=path/to/mpicc -DMPI_Fortran_COMPILER=path/to/mpif90
        # at cmake configure
        if(NOT MPI_C_COMPILER)
            set(MPI_C_COMPILER "${SIMGRID_DIR_FOUND}/bin/smpicc")
        endif()
        if(NOT MPI_C_LIBRARIES)
            set(MPI_C_LIBRARIES "${SIMGRID_LIBRARIES}")
        endif()
        if(NOT MPI_C_INCLUDE_PATH)
            set(MPI_C_INCLUDE_PATH "${SIMGRID_INCLUDE_DIRS}")
            list(APPEND MPI_C_INCLUDE_PATH "${SIMGRID_INCLUDE_DIRS}/smpi")
        endif()
#         if(NOT MPI_CXX_COMPILER)
#             set(MPI_CXX_COMPILER "${SIMGRID_DIR_FOUND}/bin/smpicxx")
#         endif()
#         if(NOT MPI_CXX_LIBRARIES)
#             set(MPI_CXX_LIBRARIES "${SIMGRID_LIBRARIES}")
#         endif()
#         if(NOT MPI_CXX_INCLUDE_PATH)
#             set(MPI_CXX_INCLUDE_PATH "${SIMGRID_INCLUDE_DIRS}")
#             list(APPEND MPI_CXX_INCLUDE_PATH "${SIMGRID_INCLUDE_DIRS}/smpi")
#         endif()
#         if(NOT MPI_Fortran_COMPILER)
#             set(MPI_Fortran_COMPILER "${SIMGRID_DIR_FOUND}/bin/smpicxx")
#         endif()
#         if(NOT MPI_Fortran_LIBRARIES)
#             set(MPI_Fortran_LIBRARIES "${SIMGRID_LIBRARIES}")
#         endif()
#         if(NOT MPI_Fortran_INCLUDE_PATH)
#             set(MPI_Fortran_INCLUDE_PATH "${SIMGRID_INCLUDE_DIRS}")
#             list(APPEND MPI_Fortran_INCLUDE_PATH "${SIMGRID_INCLUDE_DIRS}/smpi")
#         endif()
        find_package(MPI REQUIRED)

        if (MPI_C_FOUND)
            message("-- ${Blue}Add definition HICMA_USE_MPI"
            " - Activate MPI in HiCMA${ColourReset}")
            set(HICMA_USE_MPI 1)
            list(APPEND EXTRA_LIBRARIES ${MPI_C_LIBRARIES} )
            include_directories( ${MPI_C_INCLUDE_PATH} )
            if(MPI_C_LINK_FLAGS)
                list(APPEND CMAKE_EXE_LINKER_FLAGS "${MPI_C_LINK_FLAGS}")
            endif()
        endif (MPI_C_FOUND)

    endif (HICMA_USE_MPI)

endif(NOT HICMA_SIMULATION)


# HICMA depends on STARSH
# -------------------------------
find_package( STARSH REQUIRED )
if( STARSH_FOUND )
    include_directories(${STARSH_INCLUDE_DIRS_DEP})
    if(STARSH_LINKER_FLAGS)
        list(APPEND CMAKE_EXE_LINKER_FLAGS "${STARSH_LINKER_FLAGS}")
    endif()
    if(STARSH_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${STARSH_LIBRARY_DIRS}")
    endif()
    if (STARSH_LIBRARIES)
        # look for gsl
        find_library( _STARSH_LIB NAME starsh PATHS ${STARSH_LIBRARY_DIRS} )
        if( _STARSH_LIB AND NOT "${STARSH_LIBRARIES_DEP}" MATCHES "gsl" )
            execute_process( COMMAND nm ${_STARSH_LIB} COMMAND grep -q gsl RESULT_VARIABLE GSL_IN_STARSH )
            if ( ${GSL_IN_STARSH} EQUAL 0)
                message( STATUS "STARSH depends on gsl. Adding it to dependency list")
                find_package(GSL REQUIRED)
                if (STARSH_LIBRARIES_DEP)
                    list( APPEND STARSH_LIBRARIES_DEP ${GSL_LIBRARIES})
                else()
                    list( APPEND STARSH_LIBRARIES ${GSL_LIBRARIES})
                endif()
            endif()
	endif()
	# insert to dependencies
        if (STARSH_LIBRARIES_DEP)
            list(APPEND HICMA_DEP ${STARSH_LIBRARIES_DEP})
        else()
            list(APPEND HICMA_DEP ${STARSH_LIBRARIES})
        endif()
    endif()
endif()


# HICMA depends on a runtime
# -------------------------------
if( HICMA_SCHED_STARPU )

  #set(HICMA_STARPU_VERSION "1.1" CACHE STRING "oldest STARPU version desired")
  set(HICMA_STARPU_VERSION "1.2" CACHE STRING "STARPU version desired")

    # create list of components in order to make a single call to find_package(starpu...)
    if(NOT HICMA_SIMULATION)
        set(STARPU_COMPONENT_LIST "HWLOC")
#         if(HICMA_USE_CUDA)
#             list(APPEND STARPU_COMPONENT_LIST "CUDA")
#         endif()
    else()
        set(STARPU_COMPONENT_LIST "SIMGRID")
    endif()
    if(HICMA_USE_MPI)
        list(APPEND STARPU_COMPONENT_LIST "MPI")
    endif()
    if(HICMA_ENABLE_TRACING)
        list(APPEND STARPU_COMPONENT_LIST "FXT")
    endif()

    find_package(STARPU ${HICMA_STARPU_VERSION} REQUIRED
                 COMPONENTS ${STARPU_COMPONENT_LIST})

    # Add definition and include_dir if found
    if ( STARPU_FOUND )
        message("-- ${Blue}Add definition HICMA_SCHED_STARPU"
        " - Activate StarPU in HiCMA${ColourReset}")
        set(HICMA_SCHED_STARPU ON)
        include_directories(${STARPU_INCLUDE_DIRS_DEP})
        if(STARPU_LINKER_FLAGS)
            list(APPEND CMAKE_EXE_LINKER_FLAGS "${STARPU_LINKER_FLAGS}")
        endif()
        set(CMAKE_REQUIRED_INCLUDES "${STARPU_INCLUDE_DIRS_DEP}")
        foreach(libdir ${STARPU_LIBRARY_DIRS_DEP})
            list(APPEND CMAKE_REQUIRED_FLAGS "-L${libdir}")
        endforeach()
        set(CMAKE_REQUIRED_LIBRARIES "${STARPU_LIBRARIES_DEP}")
        if (HICMA_USE_MPI)
            list(APPEND CMAKE_REQUIRED_INCLUDES "${MPI_C_INCLUDE_PATH}")
            list(APPEND CMAKE_REQUIRED_FLAGS "${MPI_C_LINK_FLAGS}")
            list(APPEND CMAKE_REQUIRED_LIBRARIES "${MPI_C_LIBRARIES}")
        endif()
        if (HICMA_SIMULATION)
            list(APPEND CMAKE_REQUIRED_FLAGS "-include" "starpu_simgrid_wrap.h")
        endif()
        string(REPLACE ";" " " CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
        check_function_exists(starpu_data_idle_prefetch_on_node STARPU_IDLE_PREFETCH_FOUND)
        if ( STARPU_IDLE_PREFETCH_FOUND )
            set(HAVE_STARPU_IDLE_PREFETCH 1)
            message("-- ${Blue}Add definition HAVE_STARPU_IDLE_PREFETCH${ColourReset}")
        endif()
        check_function_exists(starpu_iteration_push STARPU_ITERATION_PUSH_FOUND)
        if ( STARPU_ITERATION_PUSH_FOUND )
            set(HAVE_STARPU_ITERATION_PUSH 1)
            message("-- ${Blue}Add definition HAVE_STARPU_ITERATION_PUSH${ColourReset}")
        endif()
        check_function_exists(starpu_data_wont_use STARPU_DATA_WONT_USE_FOUND)
        if ( STARPU_DATA_WONT_USE_FOUND )
            set(HAVE_STARPU_DATA_WONT_USE 1)
            message("-- ${Blue}Add definition HAVE_STARPU_DATA_WONT_USE${ColourReset}")
        endif()
        check_function_exists(starpu_data_set_coordinates STARPU_DATA_SET_COORDINATES_FOUND)
        if ( STARPU_DATA_SET_COORDINATES_FOUND )
            set(HAVE_STARPU_DATA_SET_COORDINATES 1)
            message("-- ${Blue}Add definition HAVE_STARPU_DATA_SET_COORDINATES${ColourReset}")
        endif()
        check_function_exists(starpu_malloc_on_node_set_default_flags STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS)
        if ( STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS )
            set(HAVE_STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS 1)
            message("-- ${Blue}Add definition HAVE_STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS${ColourReset}")
        endif()
        if(HICMA_ENABLE_TRACING)
            # check if fxt profiling is accessible in starpu and activate it in hicma
            unset(STARPU_FXT_START_PROFILING_FOUND CACHE)
            check_function_exists(starpu_fxt_start_profiling STARPU_FXT_START_PROFILING_FOUND)
            if ( STARPU_FXT_START_PROFILING_FOUND )
                message("-- ${Blue}Add definition HAVE_STARPU_FXT_PROFILING"
                " - Activate FxT profiling through StarPU${ColourReset}")
                set(HAVE_STARPU_FXT_PROFILING 1)
            else()
                message("-- ${Red}Looking for starpu with fxt"
                " - starpu_fxt_start_profiling() test fails in StarPU${ColourReset}")
                message("-- ${Red}Check in CMakeFiles/CMakeError.log to figure out why it fails${ColourReset}")
            endif()
        endif()
        if (HICMA_USE_MPI)
            # Check if a specific function exist
            unset(STARPU_MPI_DATA_REGISTER_FOUND CACHE)
            check_function_exists(starpu_mpi_data_register_comm STARPU_MPI_DATA_REGISTER_FOUND)
            if ( STARPU_MPI_DATA_REGISTER_FOUND )
                message("-- ${Blue}Add definition HAVE_STARPU_MPI_DATA_REGISTER - Activate"
                " use of starpu_mpi_data_register() in HiCMA with StarPU${ColourReset}")
                set(HAVE_STARPU_MPI_DATA_REGISTER 1)
            else()
                message("-- ${Red}Looking for starpu with starpu_mpi_data_register"
                " - starpu_mpi_data_register() test fails in StarPU${ColourReset}")
                message("-- ${Red}Check in CMakeFiles/CMakeError.log to figure out why it fails${ColourReset}")
            endif()
            unset(STARPU_MPI_COMM_RANK_FOUND CACHE)
            check_function_exists(starpu_mpi_comm_rank STARPU_MPI_COMM_RANK_FOUND)
            if ( STARPU_MPI_COMM_RANK_FOUND )
                message("-- ${Blue}Add definition HAVE_STARPU_MPI_COMM_RANK - Activate"
                " use of starpu_mpi_comm_rank() in HiCMA with StarPU${ColourReset}")
                set(HAVE_STARPU_MPI_COMM_RANK 1)
            else()
                message("-- ${Red}Looking for starpu with starpu_mpi_comm_rank"
                " - starpu_mpi_comm_rank() test fails in StarPU${ColourReset}")
                message("-- ${Red}Check in CMakeFiles/CMakeError.log to figure out why it fails${ColourReset}")
            endif()
	    check_function_exists(starpu_mpi_cached_receive STARPU_MPI_CACHED_RECEIVE)
	    if ( STARPU_MPI_CACHED_RECEIVE )
		set(HAVE_STARPU_MPI_CACHED_RECEIVE 1)
		message("-- ${Blue}Add definition HAVE_STARPU_MPI_CACHED_RECEIVE${ColourReset}")
	    endif()
        endif()
        if(HWLOC_FOUND AND HWLOC_LIBRARY_DIRS)
            # the RPATH to be used when installing
            list(APPEND CMAKE_INSTALL_RPATH "${HWLOC_LIBRARY_DIRS}")
        endif()
        if(FXT_FOUND AND FXT_LIBRARY_DIRS)
            # the RPATH to be used when installing
            list(APPEND CMAKE_INSTALL_RPATH "${FXT_LIBRARY_DIRS}")
        endif()
        if(SIMGRID_FOUND AND SIMGRID_LIBRARY_DIRS)
            # the RPATH to be used when installing
            list(APPEND CMAKE_INSTALL_RPATH "${SIMGRID_LIBRARY_DIRS}")
        endif()
        if(STARPU_FOUND AND STARPU_LIBRARY_DIRS)
            # the RPATH to be used when installing
            list(APPEND CMAKE_INSTALL_RPATH "${STARPU_LIBRARY_DIRS}")
        endif()
        if (STARPU_LIBRARIES)
            if (STARPU_LIBRARIES_DEP)
                list(INSERT HICMA_DEP 0 ${STARPU_LIBRARIES_DEP})
            else()
                list(INSERT HICMA_DEP 0 ${STARPU_LIBRARIES})
            endif()
        endif()
    else ( STARPU_FOUND )
        if(ECRC_VERBOSE_FIND_PACKAGE)
            if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
                if (NOT HWLOC_hwloc.h_DIRS)
                    Print_Find_Header_Status(hwloc hwloc.h)
                endif ()
                if (NOT HWLOC_hwloc_LIBRARY)
                    Print_Find_Library_Status(hwloc libhwloc)
                endif ()
            endif()
            if(HICMA_ENABLE_TRACING AND (NOT FXT_FOUND OR NOT FXT_LIBRARIES))
                if (NOT FXT_fxt.h_DIRS)
                    Print_Find_Header_Status(fxt fxt.h)
                endif ()
                if (NOT FXT_fxt_LIBRARY)
                    Print_Find_Library_Status(fxt libfxt)
                endif ()
            endif()
            if(HICMA_SIMULATION AND (NOT SIMGRID_FOUND OR NOT SIMGRID_LIBRARIES))
                if (NOT SIMGRID_simgrid.h_DIRS)
                    Print_Find_Header_Status(simgrid simgrid.h)
                endif ()
                if (NOT SIMGRID_simgrid_LIBRARY)
                    Print_Find_Library_Status(simgrid libsimgrid)
                endif ()
            endif()
            if( (NOT STARPU_SHM_FOUND) OR (NOT STARPU_SHM_LIBRARIES) OR
                ( STARPU_LOOK_FOR_MPI AND (NOT STARPU_MPI_FOUND OR NOT STARPU_MPI_LIBRARIES) )
              )
                foreach(starpu_hdr ${STARPU_hdrs_to_find})
                    if (NOT STARPU_${starpu_hdr}_INCLUDE_DIRS)
                        Print_Find_Header_Status(starpu ${starpu_hdr})
                    endif ()
                endforeach()
                if(STARPU_VERSION_STRING)
                    foreach(starpu_lib ${STARPU_libs_to_find})
                        if (NOT STARPU_${starpu_lib}_LIBRARY)
                            Print_Find_Library_Status(starpu ${starpu_lib})
                        endif ()
                    endforeach()
                endif ()
            endif ()
        else(ECRC_VERBOSE_FIND_PACKAGE)
            message(WARNING "StarPU library has not been found and ECRC_VERBOSE_FIND_PACKAGE is set to OFF."
            " Try to activate ECRC_VERBOSE_FIND_PACKAGE option (-DECRC_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
        endif(ECRC_VERBOSE_FIND_PACKAGE)
        if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
            message(FATAL_ERROR "hwloc library is required but has not been found")
        endif()
        if(HICMA_SIMULATION AND (NOT SIMGRID_FOUND OR NOT SIMGRID_LIBRARIES))
            message(FATAL_ERROR "SimGrid library is required but has not been found")
        endif()
        if(HICMA_ENABLE_TRACING AND (NOT FXT_FOUND OR NOT FXT_LIBRARIES))
            message(FATAL_ERROR "FxT library is required but has not been found")
        endif()
        if( (NOT STARPU_SHM_FOUND) OR (NOT STARPU_SHM_LIBRARIES) OR
            ( STARPU_LOOK_FOR_MPI AND (NOT STARPU_MPI_FOUND OR NOT STARPU_MPI_LIBRARIES) )
          )
            message(FATAL_ERROR "StarPU library is required but has not been found")
        endif()
    endif ( STARPU_FOUND )

endif( HICMA_SCHED_STARPU )

if( HICMA_SCHED_PARSEC )

    # create list of components in order to make a single call to find_package(starpu...)
    set(PARSEC_COMPONENT_LIST "HWLOC")
    if(HICMA_USE_MPI)
        list(APPEND PARSEC_COMPONENT_LIST "MPI")
    endif()
#     if(HICMA_USE_CUDA)
#         list(APPEND PARSEC_COMPONENT_LIST "CUDA")
#     endif()
    # TODO: Add a HICMA_WITH_PROFILING option that enables Fxt for StarPU, or PAPI for PaRSEC
    #if(HICMA_WITH_PROFILING)
    #    list(APPEND PARSEC_COMPONENT_LIST "PAPI")
    #endif()

    find_package(PARSEC COMPONENTS ${PARSEC_COMPONENT_LIST})

    # Add definition and include_dir if found
    if(PARSEC_FOUND)
        message("-- ${Blue}Add definition HICMA_SCHED_PARSEC"
        " - Activate PaRSEC in HiCMA${ColourReset}")

      set(HICMA_SCHED_PARSEC 1)
      include_directories(${PARSEC_INCLUDE_DIRS_DEP})
      include_directories(${PARSEC_INCLUDE_DIRS}/daguepp)
      set(CMAKE_REQUIRED_INCLUDES "${PARSEC_INCLUDE_DIRS_DEP}")

      if(PARSEC_LINKER_FLAGS)
           list(APPEND CMAKE_EXE_LINKER_FLAGS "${PARSEC_LINKER_FLAGS}")
      endif()
      foreach(libdir ${PARSEC_LIBRARY_DIRS_DEP})
          list(APPEND CMAKE_REQUIRED_FLAGS "-L${libdir}")
      endforeach()
      string(REPLACE ";" " " CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
      set(CMAKE_REQUIRED_LIBRARIES "${PARSEC_LIBRARIES_DEP}")

      if(HWLOC_FOUND AND HWLOC_LIBRARY_DIRS)
            # the RPATH to be used when installing
          list(APPEND CMAKE_INSTALL_RPATH "${HWLOC_LIBRARY_DIRS}")
      endif()
      if(PARSEC_FOUND AND PARSEC_LIBRARY_DIRS_DEP)
          # the RPATH to be used when installing
          list(APPEND CMAKE_INSTALL_RPATH "${PARSEC_LIBRARY_DIRS_DEP}")
      endif()
      if (PARSEC_LIBRARIES)
          if (PARSEC_LIBRARIES_DEP)
              list(INSERT HICMA_DEP 0 ${PARSEC_LIBRARIES_DEP})
          else()
              list(INSERT HICMA_DEP 0 ${PARSEC_LIBRARIES})
          endif()
      endif()
    else ( PARSEC_FOUND )
        if(ECRC_VERBOSE_FIND_PACKAGE)
            if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
                if (NOT HWLOC_hwloc.h_DIRS)
                    Print_Find_Header_Status(hwloc hwloc.h)
                endif ()
                if (NOT HWLOC_hwloc_LIBRARY)
                    Print_Find_Library_Status(hwloc libhwloc)
                endif ()
            endif()
            if( (NOT PARSEC_FOUND) OR (NOT PARSEC_LIBRARIES_DEP) )
              foreach(parsec_hdr ${PARSEC_hdrs_to_find})
                    if (NOT PARSEC_${parsec_hdr}_INCLUDE_DIRS)
                        Print_Find_Header_Status(parsec ${parsec_hdr})
                    endif ()
                endforeach()
                if(PARSEC_VERSION_STRING)
                    foreach(parsec_lib ${PARSEC_libs_to_find})
                        if (NOT PARSEC_${parsec_lib}_LIBRARY)
                            Print_Find_Library_Status(parsec ${parsec_lib})
                        endif ()
                    endforeach()
                endif()
            endif()
        else(ECRC_VERBOSE_FIND_PACKAGE)
            message(WARNING "PaRSEC library has not been found and ECRC_VERBOSE_FIND_PACKAGE is set to OFF."
            " Try to activate ECRC_VERBOSE_FIND_PACKAGE option (-DECRC_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
        endif(ECRC_VERBOSE_FIND_PACKAGE)
        if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
            message(FATAL_ERROR "hwloc library is required but has not been found")
        endif()
        if( (NOT PARSEC_FOUND) OR (NOT PARSEC_LIBRARIES_DEP) )
            message(FATAL_ERROR "PaRSEC library is required but has not been found")
        endif()
    endif ( PARSEC_FOUND )

endif( HICMA_SCHED_PARSEC )

# try to find quark runtime
if( HICMA_SCHED_QUARK )

    find_package(QUARK COMPONENTS HWLOC)

    # Add definition and include_dir if found
    if(QUARK_FOUND)
        message("-- ${Blue}Add definition HICMA_SCHED_QUARK"
        " - Activate QUARK in HiCMA${ColourReset}")
        set(HICMA_SCHED_QUARK 1)
        include_directories(${QUARK_INCLUDE_DIRS})
        if(QUARK_LIBRARY_DIRS)
            # the RPATH to be used when installing
            list(APPEND CMAKE_INSTALL_RPATH "${QUARK_LIBRARY_DIRS}")
        endif()
        if (QUARK_LIBRARIES)
            if (QUARK_LIBRARIES_DEP)
                list(INSERT HICMA_DEP 0 ${QUARK_LIBRARIES_DEP})
            else()
                list(INSERT HICMA_DEP 0 ${QUARK_LIBRARIES})
            endif()
        endif()
    else(QUARK_FOUND)
        if(ECRC_VERBOSE_FIND_PACKAGE)
            if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
                if (NOT HWLOC_hwloc.h_DIRS)
                    Print_Find_Header_Status(hwloc hwloc.h)
                endif ()
                if (NOT HWLOC_hwloc_LIBRARY)
                    Print_Find_Library_Status(hwloc libhwloc)
                endif ()
            endif()
            if (NOT QUARK_quark.h_DIRS)
                Print_Find_Header_Status(quark quark.h)
            endif ()
            if (NOT QUARK_quark_LIBRARY)
                Print_Find_Library_Status(quark libquark)
            endif ()
        else()
            message(WARNING "QUARK library has not been found and ECRC_VERBOSE_FIND_PACKAGE is set to OFF."
            " Try to activate ECRC_VERBOSE_FIND_PACKAGE option (-DECRC_VERBOSE_FIND_PACKAGE=ON) to get some hints for the detection")
        endif()
        if(NOT HWLOC_FOUND OR NOT HWLOC_LIBRARIES)
            message(FATAL_ERROR "hwloc library is required but has not been found")
        endif()
        message(FATAL_ERROR "QUARK library is required but has not been found")
    endif(QUARK_FOUND)

endif()

set( CHAMELEON_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/chameleon" CACHE PATH "Location of CHAMELEON source code" )
set( CHAMELEON_BUILD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/chameleon/build" CACHE PATH "Location of CHAMELEON build directory" )

#required files in source dir
include_directories( AFTER ${CHAMELEON_SOURCE_DIR}/ )


#required files in build dir
#../chameleon/bb/runtime/starpu/include/chameleon_starpu.h
include_directories( AFTER ${CHAMELEON_BUILD_DIR} )

# HICMA depends on CHAMELEON
# -------------------------------
find_package( CHAMELEON REQUIRED )
if( CHAMELEON_FOUND )
    include_directories( BEFORE ${CHAMELEON_INCLUDE_DIRS_DEP} )
    if(CHAMELEON_LINKER_FLAGS)
        list(APPEND CMAKE_EXE_LINKER_FLAGS "${CHAMELEON_LINKER_FLAGS}")
    endif()
    if(CHAMELEON_LIBRARY_DIRS)
        # the RPATH to be used when installing
        list(APPEND CMAKE_INSTALL_RPATH "${CHAMELEON_LIBRARY_DIRS}")
    endif()
    if (CHAMELEON_LIBRARIES)
        if (CHAMELEON_LIBRARIES_DEP)
            list(INSERT HICMA_DEP 0 ${CHAMELEON_LIBRARIES_DEP})
        else()
            list(INSERT HICMA_DEP 0 ${CHAMELEON_LIBRARIES})
        endif()
    endif()
endif()


# # Add option to exploit cublas API v2
# # -----------------------------------
# cmake_dependent_option(HICMA_USE_CUBLAS_V2
#   "Enable cublas API v2" ON
#   "HICMA_USE_CUDA;HICMA_SCHED_STARPU" OFF)

list(REMOVE_DUPLICATES CMAKE_EXE_LINKER_FLAGS)
string(REPLACE ";" " " CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
# Fix a problem on Mac OS X when building shared libraries
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(CMAKE_SHARED_LINKER_FLAGS "-undefined dynamic_lookup")
endif()

# # Add define for Fortran Mangling (should be defined somewhere else)
# # ------------------------------------------------------------------
# message("-- ${Blue}Add definition ADD_"
#         " - For Fortran mangling${ColourReset}")
# set(ADD_ 1)

#------------------------------------------------------------------------------
if(HICMA_SCHED_STARPU)
    link_directories(${STARPU_LIBRARY_DIRS_DEP})
elseif(HICMA_SCHED_PARSEC)
    link_directories(${PARSEC_LIBRARY_DIRS_DEP})
elseif(HICMA_SCHED_QUARK)
    link_directories(${QUARK_LIBRARY_DIRS})
endif()
if(NOT HICMA_SIMULATION)
#     if(HICMA_USE_CUDA)
#         link_directories(${CUDA_LIBRARY_DIRS})
#     endif()
#     if(HICMA_USE_MAGMA)
#         link_directories(${MAGMA_LIBRARY_DIRS})
#     endif()
endif()

# Save extra dependencies (all required links)
list(APPEND HICMA_DEP ${EXTRA_LIBRARIES})
list(REMOVE_DUPLICATES HICMA_DEP) # WARNING: is it safe, respect order?

# Check for the subdirectories
# ----------------------------

# include headers
# ---------------
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)
include_directories(${CMAKE_CURRENT_BINARY_DIR}/include)


###############################################################################
# Coreblas library (kernels for CPUs) #
#######################################

#######add_subdirectory(hcore)
#------------------------------------------------------------------------------


###############################################################################
# Cudablas library (kernels for CUDAs) #
########################################

# if(HICMA_USE_CUDA)
#     add_subdirectory(cudablas)
# endif()
#------------------------------------------------------------------------------


###############################################################################
# Main library #
################

#add_subdirectory(include)
#add_subdirectory(control)
#######add_subdirectory(compute)
#------------------------------------------------------------------------------


###############################################################################
# Runtime library #
###################

#######add_subdirectory(runtime)
#------------------------------------------------------------------------------


###############################################################################
# Executables and Tests #
#########################

enable_testing()
include(CTest)

# # Necessary to compile executables (temporary)
# # should use SystemDetection.cmake module?
# string(REGEX MATCH ".*ifort$" _match_ifort ${CMAKE_Fortran_COMPILER})
# if(_match_ifort)
#     message(STATUS "Add -nofor_main to the Fortran linker (Intel compiler)")
#     set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -nofor_main")
# endif(_match_ifort)

# Testing executables
#if(HICMA_ENABLE_EXAMPLE AND NOT HICMA_SIMULATION)
#  add_subdirectory(example)
#endif(HICMA_ENABLE_EXAMPLE AND NOT HICMA_SIMULATION)

# Testing executables
#if(HICMA_ENABLE_TESTING AND NOT HICMA_SIMULATION)
#  add_subdirectory(testing)
#endif(HICMA_ENABLE_TESTING AND NOT HICMA_SIMULATION)

# Enable MKL compilation flag
if(BLAS_VENDOR_FOUND MATCHES "Intel*")
    add_definitions(-DMKL)
endif()


# Timing executables
if(HICMA_ENABLE_TIMING)
  add_subdirectory(timing)
endif(HICMA_ENABLE_TIMING)
# Testing executables
if(HICMA_ENABLE_TESTING)
  add_subdirectory(testing)
endif(HICMA_ENABLE_TESTING)
#------------------------------------------------------------------------------

##################################################################
# temporary basic build
##################################################################
##################################################################
find_package(OpenMP )
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
else()
    ## In shaheen openmp is not detected, so i force it.
    #set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")
    #set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")
    ## It is not good to force on MAC.
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")
endif()



include_directories( ${CMAKE_SOURCE_DIR}/include )

SET( HICMA_HCORE_SRC
    ${CMAKE_SOURCE_DIR}/hcore/compute/hcore_zhagcm.c
    ${CMAKE_SOURCE_DIR}/hcore/compute/hcore_zhagdm.c
    ${CMAKE_SOURCE_DIR}/hcore/compute/hcore_zgytlr.c
    ${CMAKE_SOURCE_DIR}/hcore/compute/hcore_zsyrk.c
    ${CMAKE_SOURCE_DIR}/hcore/compute/hcore_zgemm_fast.c
    ${CMAKE_SOURCE_DIR}/hcore/compute/hcore_zgemm.c
    ${CMAKE_SOURCE_DIR}/hcore/compute/hcore_zgemmbdcd.c
    ${CMAKE_SOURCE_DIR}/hcore/compute/hcore_zuncompress.c
  )
include_directories( ${CMAKE_SOURCE_DIR}/hcore/include )

set( HICMA_CONTROL_SRC
    ${CMAKE_SOURCE_DIR}/control/hicma_init.c
    )
set( HICMA_COMPUTE_SRC
    ${CMAKE_SOURCE_DIR}/compute/pzhagcm.c
    ${CMAKE_SOURCE_DIR}/compute/zhagcm.c
    ${CMAKE_SOURCE_DIR}/compute/pzhagdm.c
    ${CMAKE_SOURCE_DIR}/compute/zhagdm.c
    ${CMAKE_SOURCE_DIR}/compute/pzgytlr.c
    ${CMAKE_SOURCE_DIR}/compute/zgytlr.c
    ${CMAKE_SOURCE_DIR}/compute/pzpotrf.c
    ${CMAKE_SOURCE_DIR}/compute/zpotrf.c
    ${CMAKE_SOURCE_DIR}/compute/pzgemm.c
    ${CMAKE_SOURCE_DIR}/compute/zgemm.c
    ${CMAKE_SOURCE_DIR}/compute/pztrsm.c
    ${CMAKE_SOURCE_DIR}/compute/ztrsm.c
  )

#set( HICMA_TIMING_SRC
    #${CMAKE_SOURCE_DIR}/timing/time_tzgemm_tile.c
    #${CMAKE_SOURCE_DIR}/timing/timing.c
  #)
include_directories( ${CMAKE_SOURCE_DIR}/timing )
include_directories( ${CMAKE_SOURCE_DIR}/testing)

set( HICMA_RUNTIME_SRC
    ${CMAKE_SOURCE_DIR}/runtime/starpu/codelets/codelet_zhagcm.c
    ${CMAKE_SOURCE_DIR}/runtime/starpu/codelets/codelet_zhagdm.c
    ${CMAKE_SOURCE_DIR}/runtime/starpu/codelets/codelet_zuncompress.c
    ${CMAKE_SOURCE_DIR}/runtime/starpu/codelets/codelet_zgemm.c
    ${CMAKE_SOURCE_DIR}/runtime/starpu/codelets/codelet_zgemm_bdcd.c
    ${CMAKE_SOURCE_DIR}/runtime/starpu/codelets/codelet_zgytlr.c
    ${CMAKE_SOURCE_DIR}/runtime/starpu/codelets/codelet_zpotrf.c
    ${CMAKE_SOURCE_DIR}/runtime/starpu/codelets/codelet_ztrsm.c
    ${CMAKE_SOURCE_DIR}/runtime/starpu/codelets/codelet_zsyrk.c
    ${CMAKE_SOURCE_DIR}/runtime/starpu/codelets/codelet_zgytlr_diag.c
  )

set( HICMA_AUX_SRC
    ${CMAKE_SOURCE_DIR}/misc/compute/zuncompress.c
    ${CMAKE_SOURCE_DIR}/misc/compute/zdiag.c
    ${CMAKE_SOURCE_DIR}/misc/compute/zproblem.c
    ${CMAKE_SOURCE_DIR}/misc/compute/znormest.c
    ${CMAKE_SOURCE_DIR}/misc/descutil/descprint.c
    ${CMAKE_SOURCE_DIR}/misc/descutil/descsame.c
    ${CMAKE_SOURCE_DIR}/misc/compute/zstat.c
  )
include_directories( ${CMAKE_SOURCE_DIR}/misc/include )

link_directories(${STARSH_LIBRARY_DIRS})
link_directories(${STARPU_LIBRARY_DIRS})
link_directories(${CHAMELEON_LIBRARY_DIRS})

add_library( hicma ${HICMA_RUNTIME_SRC} ${HICMA_HCORE_SRC} ${HICMA_COMPUTE_SRC} ${HICMA_AUX_SRC} ${HICMA_CONTROL_SRC} )
target_link_libraries(hicma ${STARSH_LIBRARIES_DEP})
#add_executable( main ${HICMA_TIMING_SRC} )
#target_link_libraries( main
    #hicma
    #${HICMA_DEP}
    #${STARSH_LIBRARIES_DEP}
    #${CHAMELEON_LIBRARIES_DEP}
    #${STARPU_LIBRARIES_DEP}
    #)

# Installation step
install( TARGETS hicma
         DESTINATION lib )
#install( TARGETS main
         #DESTINATION bin/timing )
install( FILES hcore/include/hcore_z.h
         DESTINATION include )
install( FILES include/hicma_constants.h
         DESTINATION include )
install( FILES include/hicma_struct.h
         DESTINATION include )
install( FILES include/hicma_z.h
         DESTINATION include )
install( FILES misc/include/auxcompute_z.h
         DESTINATION include/misc )
install( FILES misc/include/auxdescutil.h
         DESTINATION include/misc )

###############################################################################
# Documentation #
#################

#if(HICMA_ENABLE_DOCS)
#    add_subdirectory(docs)
#    if (DOXYGEN_EXECUTABLE AND MAKEINFO_COMPILER AND TEXI2DVI_COMPILER AND TEX_COMPILER)
#        add_custom_target(docs ALL DEPENDS doxygen-out doc-info-users_guide doc-html-users_guide doc-pdf-users_guide)
#    endif()
#endif()
#------------------------------------------------------------------------------


###############################################################################
# Config files (pkg_config) #
#############################

# Create file hicma.pc
# ---------------------------
generate_pkgconfig_file()
#------------------------------------------------------------------------------

###############################################################################
# Release source #
##################
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_NAME "hicma")
set(CPACK_PACKAGE_DESCRIPTION_FILE ${CMAKE_CURRENT_SOURCE_DIR}/README.md)
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "HiCMA is a C library providing parallel algorithms to perform BLAS/LAPACK "
                                      " operations exploiting fully modern architectures.")
set(CPACK_PACKAGE_VERSION "${HICMA_VERSION_MAJOR}.${HICMA_VERSION_MINOR}.${HICMA_VERSION_PATCH}")
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}")
set(CPACK_PACKAGE_VENDOR "Inria")
set(CPACK_PACKAGE_CONTACT "morse-devel@lists.gforge.inria.fr")
set(CPACK_RESOURCE_FILE_README ${CMAKE_CURRENT_SOURCE_DIR}/README.md)
set(CPACK_RESOURCE_FILE_LICENSE ${CMAKE_CURRENT_SOURCE_DIR}/LICENSE)
set(CPACK_SOURCE_IGNORE_FILES "plasma-conversion;build;.cproject;.settings;.dir-locals.el;.project;.pydevproject;.svn;.git;.gitmodules;.gitlab-ci.yml")
include(CPack)

###############################################################################
# Print Options #
#################
get_directory_property( HICMA_DEFINITIONS_LIST DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS )
include(PrintOpts)

if (CMAKE_INSTALL_PREFIX STREQUAL "/usr/local" OR CMAKE_INSTALL_PREFIX STREQUAL "C:/Program Files")
    message("-- ${Yellow}Your CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX} which is a default system path."
    " You may want to change it: set the CMAKE_INSTALL_PREFIX variable to do so${ColourReset}")
else()
    message("-- ${Yellow}Your CMAKE_INSTALL_PREFIX is ${CMAKE_INSTALL_PREFIX}${ColourReset}")
endif()
message("VERY IMPORTANT")
message("Current build system assumes that Chameleon is built at ${CHAMELEON_BUILD_DIR}")
message("This cmake file can know install dir: ${CHAMELEON_DIR_FOUND}")
message("This cmake file can NOT know build dir!!!!!!!!!!")

###
### END CMakeLists.txt
###
