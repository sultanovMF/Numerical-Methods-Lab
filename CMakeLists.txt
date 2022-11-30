cmake_minimum_required (VERSION 3.8)

project ("Numerical Analysis")

set(BUILD_SHARED_LIBS OFF)

find_program(CCACHE_FOUND ccache)
if(CCACHE_FOUND)
  set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE ccache)
endif(CCACHE_FOUND)

add_subdirectory(labs)
add_subdirectory(murlib)
add_subdirectory(dependencies)