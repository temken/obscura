add_executable(tests
    test_Astronomy.cpp
    test_Configuration.cpp
    test_Direct_Detection.cpp
    test_Direct_Detection_Ionization.cpp
    test_Direct_Detection_Nucleus.cpp
    test_Direct_Detection_Semiconductor.cpp
    test_DM_Distribution.cpp
    test_DM_Particle.cpp
    test_DM_Particle_Standard.cpp
    test_Target_Electron.cpp
    test_Target_Nucleus.cpp
    )

# Linking up all libraries (includes the coverage settings set by 'example' library)
target_link_libraries(tests 
	PRIVATE
		libobscura
		gtest_main	#contains the main function
)

target_compile_options(tests PUBLIC -Wall -pedantic)

target_include_directories(tests
    PRIVATE
        ${GENERATED_DIR} )

install(TARGETS tests DESTINATION ${TESTS_DIR})

add_test(NAME obscura-tests COMMAND tests
			WORKING_DIRECTORY ${TESTS_DIR})