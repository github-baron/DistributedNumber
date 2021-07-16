cmake_minimum_required(VERSION 3.10)

#set the project name
project(DistributedNumber_Test)

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

#add the executable
#http://www.cmake.org/pipermail/cmake/2006-October/011446.html
#Modified to use pkg config and use standard var names
add_executable(DistributedNumber_Test DistributedNumber_TestRunner.cpp)


#Find the CppUnit includes and library

INCLUDE(FindPkgConfig)
PKG_CHECK_MODULES(PC_CPPUNIT "cppunit")

FIND_PATH(FLOATING_MEASURE_INCLUDE_DIRS
    NAMES cppunit/TestCase.h
    HINTS ${PROJECT_SOURCE_DIR}/../FloatingMeasure
    ${CMAKE_INSTALL_PREFIX}/../build/CMAKE/include 
)

FIND_PATH(CPPUNIT_INCLUDE_DIRS
    NAMES cppunit/TestCase.h
    HINTS ${PC_CPPUNIT_INCLUDE_DIR}
    ${CMAKE_INSTALL_PREFIX}/include
    PATHS
    /usr/local/include
    /usr/include/
    C:\\Users\\Michae\\Source\\Repos\\Ultimaker\\CppUnit\\src\\cppunit\\include
)


message("found cppunit includes " ${CPPUNIT_INCLUDE_DIRS}) 
include_directories(../DistributedNumber/ ${CPPUNIT_INCLUDE_DIRS})

FIND_LIBRARY(CPPUNIT_LIBRARIES
    NAMES cppunit
    HINTS ${PC_CPPUNIT_LIBDIR}
#    ${CMAKE_INSTALL_PREFIX}/lib
#    ${CMAKE_INSTALL_PREFIX}/lib64
#    PATHS
#    ${CPPUNIT_INCLUDE_DIRS}/../lib
    /usr/local/lib
    /usr/lib
    /usr/lib/x86_64-linux-gnu/
    C:\\Users\\Michae\\Source\\Repos\\Ultimaker\\CppUnit\\src\\cppunit\\lib
)

message("found cppunit libraries " ${CPPUNIT_LIBRARIES}) 

INCLUDE(FindPackageHandleStandardArgs)
#FIND_PACKAGE_HANDLE_STANDARD_ARGS(CPPUNIT DEFAULT_MSG CPPUNIT_LIBRARIES CPPUNIT_INCLUDE_DIRS)
MARK_AS_ADVANCED(CPPUNIT_LIBRARIES CPPUNIT_INCLUDE_DIRS)


#set version to be useable by calling script

IF(CPPUNIT_FOUND)
    set(CPPUNIT_VERSION ${PC_CPPUNIT_VERSION} CACHE INTERNAL "CppUnit Version" FORCE)
ENDIF()


# link cppunit and DistributedNumber lib
target_link_libraries(DistributedNumber_Test ${CPPUNIT_LIBRARIES} ${LOG4CPLUS_LIBRARIES} DistributedNumber)

#install into cmake / bin
# installation directory
set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/../build/CMake/)
install(TARGETS DistributedNumber_Test RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
