add_executable(tests tests.cpp)

# Linking up all libraries (includes the coverage settings set by 'example' library)
target_link_libraries(tests 
						PRIVATE obscura
						PRIVATE libphys)

add_test(NAME obscura-tests COMMAND tests
			WORKING_DIRECTORY ${TESTS_DIR})

install(TARGETS tests DESTINATION ${TESTS_DIR})
