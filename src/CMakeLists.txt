#Tester program for the new implementation
set(DWELL_REGION_EXISTS_SOURCES dwell_region_exists.cpp)
add_executable(dwell_region_exists "${DWELL_REGION_EXISTS_SOURCES}")

IF(UNIX) 
  target_link_libraries(dwell_region_exists PRIVATE gsl-lite dwell_regions ${Boost_LIBRARIES})
ENDIF(UNIX)

IF(WIN32) 
  target_link_libraries(dwell_region_exists PRIVATE gsl-lite dwell_regions ${Boost_LIBRARIES})
ENDIF(WIN32)

#Tester program for the offline dwell region queries
set(OFFLINE_DWELL_REGION_EXISTS_SOURCE offline_dwell_region_exists.cpp)
add_executable(offline_dwell_region_exists "${OFFLINE_DWELL_REGION_EXISTS_SOURCE}")

IF(UNIX) 
  target_link_libraries(offline_dwell_region_exists PRIVATE gsl-lite dwell_regions ${Boost_LIBRARIES}) 
ENDIF(UNIX)

IF(WIN32) 
  target_link_libraries(offline_dwell_region_exists PRIVATE gsl-lite dwell_regions ${Boost_LIBRARIES}) 
ENDIF(WIN32)
