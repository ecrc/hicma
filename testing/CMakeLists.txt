
set( TESTINGS_SRC 
    testing_zposv.c 
    testing_ztrsmd.c 
   )


link_directories(${STARSH_LIBRARY_DIRS})
link_directories(${STARPU_LIBRARY_DIRS})
link_directories(${CHAMELEON_LIBRARY_DIRS})



#message(${CHAMELEON_LIBRARY_DIRS})
#link_directories(${HICMA_LIBRARY_DIRS})
set(TESTING_AUX_HDRS
    ../timing/timing.h
    testing_zauxiliary.c
    testing_zauxiliary.h
    )

set(libs_for_testings)
list(APPEND libs_for_testings hicma chameleon)

add_custom_target(testing_include ALL SOURCES ${TESTING_AUX_HDRS})
include_directories( ${CMAKE_SOURCE_DIR}/misc/include )
include_directories( ${CMAKE_SOURCE_DIR}/chameleon/ )
include_directories(BEFORE ${CMAKE_SOURCE_DIR}/timing )
foreach(_timing ${TESTINGS_SRC})
    get_filename_component(_name_exe ${_timing} NAME_WE)
    add_executable(${_name_exe} ${_timing} ${TESTING_AUX_HDRS})
    add_dependencies(${_name_exe} testing_include)
    #    set_property(TARGET ${_name_exe} PROPERTY LINKER_LANGUAGE Fortran)
    target_link_libraries(${_name_exe} ${libs_for_testings})
    target_link_libraries(${_name_exe}
        hicma
        ${HICMA_DEP}
        ${STARSH_LIBRARIES_DEP}
        ${CHAMELEON_LIBRARIES_DEP}
        ${STARPU_LIBRARIES_DEP}
        )
    install(TARGETS ${_name_exe}
        DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/hicma/testing)

    # Test
    if( _name_exe MATCHES "zposv")
        add_test( NAME ${_name_exe}
                  COMMAND ${_name_exe} 1 0 posv 1000 1000 200 1000 250 1e-7 0 250 250 1 1 1)
        set_tests_properties( ${_name_exe} PROPERTIES ENVIRONMENT "MORSE_TESTING_VERBOSE=1" )
        set_tests_properties( ${_name_exe} PROPERTIES LABELS "testing" )
    endif()

endforeach()

