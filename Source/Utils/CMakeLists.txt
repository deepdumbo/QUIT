option( BUILD_UTILITIES "Build the utility programs" ON )
if( ${BUILD_UTILITIES} )
    set( PROGRAMS
        qicomplex qiaffine qimask qikfilter
        qisplitsubjects qipolyfit qipolyimg
        qi_coil_combine qi_rfprofile qi_gradient )

    foreach(PROGRAM ${PROGRAMS})
        add_executable(${PROGRAM} ${PROGRAM}.cpp)
        target_link_libraries(${PROGRAM} qi_imageio qi_core ${ITK_LIBRARIES} ${CERES_LIBRARIES})
    endforeach(PROGRAM)

    install( TARGETS ${PROGRAMS} RUNTIME DESTINATION bin )
endif()