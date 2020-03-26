file(GLOB SOURCES ${SRC_DIR}/*.cpp)

add_executable(DM_Detection ${SOURCES})

target_link_libraries(DM_Detection
	PUBLIC coverage_config 
	PUBLIC libphys)

add_test(NAME DM_Detection_test COMMAND DM_Detection)

install(TARGETS DM_Detection DESTINATION ${BIN_DIR})

