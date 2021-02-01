enable_testing()

# Testing the heap implementation
set(HEAP_TEST_SOURCES test_heap.cpp)
add_executable(heap_test ${HEAP_TEST_SOURCES} )
target_link_libraries(heap_test dwell_regions)
add_test(heap_test heap_test)

add_executable(welzl_test "test_welzl.cpp")
target_link_libraries(welzl_test dwell_regions)
add_test(welzl_test welzl_test)

add_executable(disk_intersection_test "test_disk_intersection.cpp")
target_link_libraries(disk_intersection_test dwell_regions)
add_test(disk_intersection_test disk_intersection_test)
