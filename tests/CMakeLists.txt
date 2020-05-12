add_executable(tests
	test_Astronomy.cpp)

# Linking up all libraries (includes the coverage settings set by 'example' library)
target_link_libraries(tests 
	PRIVATE
		libobscura
		gtest_main	#contains the main function
)

target_compile_options(tests PUBLIC -Wall -pedantic)

install(TARGETS tests DESTINATION ${TESTS_DIR})

add_test(NAME obscura-tests COMMAND tests
			WORKING_DIRECTORY ${TESTS_DIR})