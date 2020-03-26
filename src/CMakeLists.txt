file(GLOB SOURCES ${SRC_DIR}/*.cpp)

add_executable(executable ${SOURCES})

target_link_libraries(executable
	PUBLIC coverage_config 
	PUBLIC libphys)

add_test(NAME executable_test COMMAND executable)

install(TARGETS executable DESTINATION ${BIN_DIR})

