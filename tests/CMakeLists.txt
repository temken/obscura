add_executable(tests tests.cpp)

# Linking up all libraries (includes the coverage settings set by 'example' library)
target_link_libraries(tests 
						PRIVATE darkmatter
						PRIVATE libphys)

add_test(NAME DM_Detection_test COMMAND tests
			WORKING_DIRECTORY ${TESTS_DIR})

install(TARGETS tests DESTINATION ${TESTS_DIR})
