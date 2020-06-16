
set( TIMINGS_SRC 
    time_zgemm_tile.c 
    time_zpotrf_tile.c 
    time_zpotrf_tile_batch.c
   )


link_directories(${STARSH_LIBRARY_DIRS})
link_directories(${STARPU_LIBRARY_DIRS})
link_directories(${CHAMELEON_LIBRARY_DIRS})



#message(${CHAMELEON_LIBRARY_DIRS})
#link_directories(${HICMA_LIBRARY_DIRS})
set(TIMING_AUX_HDRS
    timing.h
    timing.c
    timing_zauxiliary.h
    timing_zauxiliary.c
   )

set(libs_for_timings)
list(APPEND libs_for_timings hicma chameleon)

add_custom_target(timing_include ALL SOURCES ${TIMING_AUX_HDRS})
include_directories( ${CMAKE_SOURCE_DIR}/misc/include )
include_directories( ${CMAKE_SOURCE_DIR}/chameleon/ )
foreach(_timing ${TIMINGS_SRC})
    get_filename_component(_name_exe ${_timing} NAME_WE)
    add_executable(${_name_exe} ${_timing} ${TIMING_AUX_HDRS})
    add_dependencies(${_name_exe} timing_include)
    #    set_property(TARGET ${_name_exe} PROPERTY LINKER_LANGUAGE Fortran)
    target_link_libraries(${_name_exe} ${libs_for_timings})
    target_link_libraries(${_name_exe}
        hicma
        ${HICMA_DEP}
        ${STARSH_LIBRARIES_DEP}
        ${CHAMELEON_LIBRARIES_DEP}
        ${STARPU_LIBRARIES_DEP}
        )
    install(TARGETS ${_name_exe}
        DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/hicma/timing)

    # Test
    if( _name_exe STREQUAL "time_zpotrf_tile")
        if( MPI_FOUND )
            add_test( NAME ${_name_exe}-ss
                      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${_name_exe} --m=3600 --n_range=3600:3600 --k=3600 --mb=400 --nb=400 --nowarmup --threads=4 --rk=0 --acc=1e-8 --check --ss --starshdecay=2 --starshmaxrank=400 )
            set_tests_properties( ${_name_exe}-ss PROPERTIES LABELS "mpi;timing" )

            add_test( NAME ${_name_exe}-edsin
                      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${_name_exe} --m=3600 --n_range=3600:3600 --k=3600 --mb=400 --nb=400 --nowarmup --threads=4 --rk=0 --acc=1e-8 --check --edsin --starshwavek=40 --starshdecay=2 --starshmaxrank=400 )
            set_tests_properties( ${_name_exe}-edsin PROPERTIES LABELS "mpi;timing" )
            
            add_test( NAME ${_name_exe}-mpi-st-3d-exp
                      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${_name_exe} --m=2500 --n_range=2500:2500 --k=2500 --mb=500 --nb=500 --nowarmup --threads=4 --check --starshmaxrank=500 --st-3D-exp --acc=1e-5 )
            set_tests_properties( ${_name_exe}-mpi-st-3d-exp PROPERTIES LABELS "mpi;timing" )
            set_tests_properties( ${_name_exe}-mpi-st-3d-exp PROPERTIES ENVIRONMENT "STARPU_SILENT=1" )
            add_test( NAME ${_name_exe}-m-3D-rbf1
                      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${_name_exe} --m=10370 --n_range=10370:10370 --k=10370 --mb=1037 --nb=1037 --nowarmup --threads=8 --rk=0 --acc=1e-6 --check --m-3D-rbf --starshdecay=0 --starshmaxrank=2047 --starshwavek=0 --rbf_kernel=9 --rad=0.6 --numobj=1 --order=2 --mesh_file=${CMAKE_SOURCE_DIR}/stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus10370.txt  --p=1)
            set_tests_properties( ${_name_exe}-m-3D-rbf1 PROPERTIES LABELS "mpi;timing" )
            set_tests_properties( ${_name_exe}-m-3D-rbf1 PROPERTIES ENVIRONMENT "STARPU_SILENT=1" )
           # add_test( NAME ${_name_exe}-m-3D-rbf2
           #           COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${_name_exe} --m=103700 --n_range=3000:3000 --k=103700 --mb=2074 --nb=60 --nowarmup --threads=8 --rk=0 --acc=1e-6  --m-3D-rbf --starshdecay=0 --starshmaxrank=120 --starshwavek=0 --rbf_kernel=9 --rad=-1 --numobj=10  --order=2 --mesh_file=${CMAKE_SOURCE_DIR}/stars-h/SARS-CoV-2-meshes/S1data/SortVirus103700.txt  --p=1)
            #set_tests_properties( ${_name_exe}-m-3D-rbf2 PROPERTIES LABELS "mpi;timing" )
            #set_tests_properties( ${_name_exe}-m-3D-rbf2 PROPERTIES ENVIRONMENT "STARPU_SILENT=1" )

#            add_test( NAME ${_name_exe}-m-3D-rbf3
#                      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${_name_exe} --m=207400 --n_range=6800:6800 --k=207400 --mb=3050 --nb=100 --nowarmup --threads=8 --rk=0 --acc=1e-6  --m-3D-rbf --starshdecay=0 --starshmaxrank=200 --starshwavek=0 --rbf_kernel=9 --rad=-1 --numobj=10  --order=2 --mesh_file=${CMAKE_SOURCE_DIR}/stars-h/SARS-CoV-2-meshes/S1data/SortVirus207400.txt  --p=2)
#            set_tests_properties( ${_name_exe}-m-3D-rbf3 PROPERTIES LABELS "mpi;timing" )
#            set_tests_properties( ${_name_exe}-m-3D-rbf3 PROPERTIES ENVIRONMENT "STARPU_SILENT=1" )
        else()
            add_test( NAME ${_name_exe}-ss
                      COMMAND ${_name_exe} --m=3600 --n_range=3600:3600 --k=3600 --mb=400 --nb=400 --nowarmup --threads=4 --rk=0 --acc=1e-8 --check --ss --starshdecay=2 --starshmaxrank=400 )
            set_tests_properties( ${_name_exe}-ss PROPERTIES LABELS "timing" )

            add_test( NAME ${_name_exe}-edsin
                      COMMAND ${_name_exe} --m=3600 --n_range=3600:3600 --k=3600 --mb=400 --nb=400 --nowarmup --threads=4 --rk=0 --acc=1e-8 --check --edsin --starshdecay=2 --starshmaxrank=400 --starshwavek=40 )
            set_tests_properties( ${_name_exe}-edsin PROPERTIES LABELS "timing" )
            add_test( NAME ${_name_exe}-m-3D-rbf1
                      COMMAND ${_name_exe} --m=10370 --n_range=10370:10370 --k=10370 --mb=1037 --nb=1037 --nowarmup --threads=16 --rk=0 --acc=1e-6 --check --m-3D-rbf --starshdecay=0 --starshmaxrank=2047 --starshwavek=0 --rbf_kernel=9 --rad=0.6 --numobj=1 --order=2 --mesh_file=${CMAKE_SOURCE_DIR}/stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus10370.txt  --p=1)
            set_tests_properties( ${_name_exe}-m-3D-rbf1 PROPERTIES LABELS "timing" )
           # add_test( NAME ${_name_exe}-m-3D-rbf2
           #           COMMAND ${_name_exe} --m=103700 --n_range=3000:3000 --k=103700 --mb=2074 --nb=60 --nowarmup --threads=16 --rk=0 --acc=1e-6  --m-3D-rbf --starshdecay=0 --starshmaxrank=120 --starshwavek=0 --rbf_kernel=9 --rad=-1 --numobj=10 --order=2 --mesh_file=${CMAKE_SOURCE_DIR}/stars-h/SARS-CoV-2-meshes/S1data/SortVirus103700.txt  --p=1)
           # set_tests_properties( ${_name_exe}-m-3D-rbf2 PROPERTIES LABELS "timing" )
           # add_test( NAME ${_name_exe}-m-3D-rbf3
           #           COMMAND ${_name_exe} --m=117715 --n_range=3900:3900 --k=117715 --mb=1811 --nb=60 --nowarmup --threads=16 --rk=0 --acc=1e-6  --m-3D-rbf --starshdecay=0 --starshmaxrank=120 --starshwavek=0 --rbf_kernel=9 --rad=0.6 --numobj=1 --order=2 --mesh_file=${CMAKE_SOURCE_DIR}/stars-h/SARS-CoV-2-meshes/singleviursdata/SortVirus117715.txt  --p=1)
           # set_tests_properties( ${_name_exe}-m-3D-rbf3 PROPERTIES LABELS "timing" )
            #add_test( NAME ${_name_exe}-m-3D-rbf4
            #          COMMAND ${_name_exe} --m=207400 --n_range=6800:6800 --k=207400 --mb=3050 --nb=100 --nowarmup --threads=16 --rk=0 --acc=1e-6  --m-3D-rbf --starshdecay=0 --starshmaxrank=200 --starshwavek=0 --rbf_kernel=9 --rad=-1 --numobj=10 --order=2 --mesh_file=${CMAKE_SOURCE_DIR}/stars-h/SARS-CoV-2-meshes/S1data/SortVirus207400.txt  --p=1)
            #set_tests_properties( ${_name_exe}-m-3D-rbf4 PROPERTIES LABELS "timing" )
        endif()
        set_tests_properties( ${_name_exe}-ss PROPERTIES ENVIRONMENT "STARPU_SILENT=1" )
        set_tests_properties( ${_name_exe}-edsin PROPERTIES ENVIRONMENT "STARPU_SILENT=1" )
    endif()
    if( _name_exe STREQUAL "time_zgemm_tile")
        if( MPI_FOUND )
            add_test( NAME ${_name_exe}-mpi-rnd
                      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${_name_exe} --m=360 --n_range=360:360 --k=360 --mb=90 --nb=90 --nowarmup --threads=1 --rk=0 --acc=1e-3 --check --rnd --starshdecay=2 --starshmaxrank=90 )
            set_tests_properties( ${_name_exe}-mpi-rnd PROPERTIES LABELS "mpi;timing" )
            add_test( NAME ${_name_exe}-mpi-edsin
                      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${_name_exe} --p=2 --m=3600 --n_range=3600:3600 --k=3600 --mb=400 --nb=400 --nowarmup --threads=10 --rk=0 --acc=1e-3 --check --edsin --starshdecay=2 --starshmaxrank=300 --starshwavek=1 )
            set_tests_properties( ${_name_exe}-mpi-edsin PROPERTIES LABELS "mpi;timing" )
            add_test( NAME ${_name_exe}-mpi-ss
                      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${_name_exe} --m=400 --n_range=400:400 --k=400 --mb=100 --nb=100 --nowarmup --threads=3 --check --starshmaxrank=100 --ss --acc=1e-2 )
            set_tests_properties( ${_name_exe}-mpi-ss PROPERTIES LABELS "mpi;timing" )
            #add_test( NAME ${_name_exe}-mpi-st-3d-exp
            #          COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${_name_exe} --m=2500 --n_range=2500:2500 --k=2500 --mb=500 --nb=500 --nowarmup --threads=3 --check --starshmaxrank=500 --st-3D-exp --acc=1e-5 )
            #set_tests_properties( ${_name_exe}-mpi-st-3d-exp PROPERTIES LABELS "mpi;timing" )

            set_tests_properties( ${_name_exe}-mpi-rnd PROPERTIES ENVIRONMENT "STARPU_SILENT=1" )
            set_tests_properties( ${_name_exe}-mpi-edsin PROPERTIES ENVIRONMENT "STARPU_SILENT=1" )
            set_tests_properties( ${_name_exe}-mpi-ss PROPERTIES ENVIRONMENT "STARPU_SILENT=1" )
            #set_tests_properties( ${_name_exe}-mpi-st-3d-exp PROPERTIES ENVIRONMENT "STARPU_SILENT=1" )
        else()
            add_test( NAME ${_name_exe}-rnd
                      COMMAND ${_name_exe} --m=360 --n_range=360:360 --k=360 --mb=90 --nb=90 --nowarmup --threads=1 --rk=0 --acc=1e-3 --check --rnd --starshdecay=2 --starshmaxrank=90 )
            set_tests_properties( ${_name_exe}-rnd PROPERTIES LABELS "timing" )
            add_test( NAME ${_name_exe}-edsin
                      COMMAND ${_name_exe} --m=3600 --n_range=3600:3600 --k=3600 --mb=400 --nb=400 --nowarmup --threads=10 --rk=0 --acc=1e-3 --check --edsin --starshdecay=2 --starshmaxrank=300 --starshwavek=1 )
            set_tests_properties( ${_name_exe}-edsin PROPERTIES LABELS "mpi;timing" )
            
            set_tests_properties( ${_name_exe}-rnd PROPERTIES ENVIRONMENT "STARPU_SILENT=1" )
            set_tests_properties( ${_name_exe}-edsin PROPERTIES ENVIRONMENT "STARPU_SILENT=1" )
        endif()
    endif()

    if( _name_exe STREQUAL "time_zpotrf_tile_batch")
        if( MPI_FOUND )
            add_test( NAME ${_name_exe}-m-3D-rbf5
                      COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${_name_exe} --m=20740 --n_range=1000:1000 --k=20740 --mb=1037 --nb=50 --nowarmup --threads=8 --rk=0 --acc=1e-6  --m-3D-rbf --starshdecay=0 --starshmaxrank=100 --starshwavek=0 --rbf_kernel=9 --rad=-1  --numobj=8  --numsubobj=2 --order=2 --mesh_file=${CMAKE_SOURCE_DIR}/stars-h/SARS-CoV-2-meshes/S2data_20k/  --p=1)
            set_tests_properties( ${_name_exe}-m-3D-rbf5 PROPERTIES LABELS "mpi;timing" )
            set_tests_properties( ${_name_exe}-m-3D-rbf5 PROPERTIES ENVIRONMENT "STARPU_SILENT=1" )
        endif()
    endif()

endforeach()
