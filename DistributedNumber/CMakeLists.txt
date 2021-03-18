cmake_minimum_required(VERSION 3.0.0)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)

# setting variables: 
# need compile commands for clazy analyzing the code
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# the source files that will create the library,
set(distributednumber_lib
    CartesianCoordinate.cpp
    Distribution.cpp
    ProbabilityDensityDistribution.cpp
)
add_library(DistributedNumber SHARED ${distributednumber_lib})

if(EXISTS "${CMAKE_SOURCE_DIR}/.git")
	execute_process(
		COMMAND git rev-parse HEAD
		WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
		OUTPUT_VARIABLE GIT_BRANCH
		OUTPUT_STRIP_TRAILING_WHITESPACE
	)
else()
	set(GIT_BRANCH "Unknown")
endif()
add_definitions(-DGIT_BRANCH=${GIT_BRANCH})
message( STATUS "LIBRARY_OUTPUT_PATH: " ${LIBRARY_OUTPUT_PATH} )



# installation directory
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/build/CMake/)
install(TARGETS DistributedNumber DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(DIRECTORY ${PROJECT_SOURCE_DIR}/DistributedNumber/ DESTINATION ${CMAKE_INSTALL_PREFIX}/include
        FILES_MATCHING PATTERN "*.h")
