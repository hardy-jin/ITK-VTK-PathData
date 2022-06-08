cmake_minimum_required(VERSION 2.8)

project (VTK_asn)

find_package (ITK REQUIRED)
include (${ITK_USE_FILE})

find_package (VTK REQUIRED)
include (${VTK_USE_FILE})

add_executable (VTK_asn VTK_asn.cxx)

target_link_libraries (VTK_asn ${ITK_LIBRARIES} ${VTK_LIBRARIES})