cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(DwellRegions LANGUAGES CXX)

cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)
if (_host_name STREQUAL breeze)
	set(CMAKE_CXX_COMPILER "/usr/local/bin/g++")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} I/usr/local/bin/include")
endif()

set(CMAKE_CXX_STANDARD 20)

IF (UNIX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -pedantic -fdiagnostics-color=always")
	set(CMAKE_CXX_FLAGS_DEBUG "-g")
	set(CMAKE_CXX_FLAGS_RELEASE "-O3 -g")
ENDIF(UNIX)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/")

# Add dependencies
# 1. GeographicLib
IF (UNIX)
	find_package (GeographicLib REQUIRED)
ELSEIF (WIN32)
	find_package (GeographicLib CONFIG REQUIRED)
ENDIF()
include_directories (${GeographicLib_INCLUDE_DIRS})
add_definitions (${GeographicLib_DEFINITIONS})

# 2. Gsl-Lite
add_library( gsl-lite INTERFACE )
target_include_directories(gsl-lite INTERFACE include)

# 3. Boost Program options
FIND_PACKAGE( Boost 1.66	COMPONENTS program_options REQUIRED )
INCLUDE_DIRECTORIES( ${Boost_INCLUDE_DIR} )

# Build code
include_directories("include/")

add_subdirectory("lib/")
add_subdirectory("src/")
add_subdirectory("tests/")
#add_subdirectory("server/")
