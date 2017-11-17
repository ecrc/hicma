###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2016 Inria. All rights reserved.
# @copyright (c) 2012-2014, 2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
# @copyright (c) 2017, King Abdullah University of Science and Technology. All rights reserved.
#
###
#
#  @file PrintOpts.cmake
#
#  @project MORSE
#  MORSE is a software package provided by:
#     Inria Bordeaux - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @author Florent Pruvost
#  @author Eduardo Gonzalez
#  @date 22-08-2017
#
###
set(dep_message "\nConfiguration of HiCMA:\n"
        "       Compiler: C .........: ${CMAKE_C_COMPILER} (${CMAKE_C_COMPILER_ID})\n"
#        "       Compiler: Fortran ...: ${CMAKE_Fortran_COMPILER} (${CMAKE_Fortran_COMPILER_ID})\n"
)
if(HICMA_USE_MPI)
  set(dep_message "${dep_message}"
  "       Compiler: MPI .......: ${MPI_C_COMPILER}\n"
  "       compiler flags ......: ${MPI_C_COMPILE_FLAGS}\n")
endif()
set(dep_message "${dep_message}"
"       Linker: .............: ${CMAKE_LINKER}\n"
"\n"
"       Build type ..........: ${CMAKE_BUILD_TYPE}\n"
"       Build shared ........: ${BUILD_SHARED_LIBS}\n"
"       CFlags ..............: ${CMAKE_C_FLAGS}\n"
"       LDFlags (shared lib).: ${CMAKE_SHARED_LINKER_FLAGS}\n"
"       LDFlags (static lib).: ${CMAKE_STATIC_LINKER_FLAGS}\n"
"       LDFlags (executable).: ${CMAKE_EXE_LINKER_FLAGS}\n"
"\n"
"       Implementation paradigm\n"
"       CUDA ................: ${HICMA_USE_CUDA}\n"
"       MPI .................: ${HICMA_USE_MPI}\n"
"\n"
"       Runtime specific\n"
"       QUARK ...............: ${HICMA_SCHED_QUARK} ${QUARK_DIR_FOUND}\n"
"       StarPU ..............: ${HICMA_SCHED_STARPU} ${STARPU_DIR_FOUND}\n"
"\n"
"       Kernels specific\n"
"       BLAS ................: ${BLAS_VENDOR_FOUND} [${BLAS_LIBRARIES}]\n"
"       LAPACK...............: ${LAPACK_VENDOR_FOUND} [${LAPACK_LIBRARIES}]\n"
"\n"
"       Chameleon ...........: ${CHAMELEON_DIR_FOUND}\n"
"       STARS-H .............: ${STARSH_DIR_FOUND}\n"
"\n"
"       Trace ...............: ${HICMA_ENABLE_TRACING}\n"
"       Simulation mode .....: ${HICMA_SIMULATION}\n"
"\n"
"       Binaries to build\n"
"       documentation ........: ${HICMA_ENABLE_DOCS}\n"
"       example ..............: ${HICMA_ENABLE_EXAMPLE}\n"
"       testing ..............: ${HICMA_ENABLE_TESTING}\n"
"       timing ...............: ${HICMA_ENABLE_TIMING}\n"
"\n"
"       HICMA dependencies :\n")
foreach (_dep ${HICMA_DEP})
    set(dep_message "${dep_message}"
    "                                 ${_dep}\n")
endforeach ()
string(REGEX REPLACE ";" " " HICMA_DEFINITIONS_LIST "${HICMA_DEFINITIONS_LIST}")
set(dep_message "${dep_message}"
"\n"
"       Definitions: ${HICMA_DEFINITIONS_LIST}\n")
set(dep_message "${dep_message}"
"\n"
"       INSTALL_PREFIX ......: ${CMAKE_INSTALL_PREFIX}\n\n")

string(REPLACE ";" " " dep_message_wsc "${dep_message}")
message(${dep_message})
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/config.log "${dep_message_wsc}")
message(STATUS "Configuration is done - A summary of the current configuration"
"\n   has been written in ${CMAKE_CURRENT_BINARY_DIR}/config.log")
# installation
# ------------
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/config.log DESTINATION share/hicma)
