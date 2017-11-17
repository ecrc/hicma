
set( TIMINGS_SRC 
    time_zgemm_tile.c 
    time_zpotrf_tile.c 
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
include_directories( ${CMAKE_SOURCE_DIR}/aux/include )
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
endforeach()

