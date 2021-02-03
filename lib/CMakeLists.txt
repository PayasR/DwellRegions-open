#Build the library containing sources for dwell region computation

set(DWELL_REGION_LIBRARY_SOURCES dwell_region_computation_base.cpp
								online_dwell_region_computation.cpp
								offline_dwell_region_rho_index.cpp
								offline_dwell_region_query_rho_index.cpp
								offline_dwell_region_tau_index.cpp
								offline_dwell_region_query_tau_index.cpp
								offline_dwell_region_query_no_index.cpp
								point.cpp
								disk.cpp
								util.cpp)
add_library(dwell_regions STATIC "${DWELL_REGION_LIBRARY_SOURCES}")
if (UNIX)
	target_link_libraries(dwell_regions ${GeographicLib_LIBRARIES})
endif(UNIX)

if (WIN32)
	target_link_libraries(dwell_regions ${GeographicLib_LIBRARIES} GeographicLib)
endif(WIN32)
