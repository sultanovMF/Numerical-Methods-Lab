project (murlib)

add_library(${PROJECT_NAME} 
"src/miscellaneous.cpp" 
"include/murlib/matrix_solver.h" 
"include/murlib/approximation.h" 
"include/murlib/integration.h" 
"src/interpolation.cpp" "src/matrix_solver.cpp" "src/approximation.cpp" "src/integration.cpp" "include/murlib/miscellaneous.h")
add_library(murlib::interpolation ALIAS ${PROJECT_NAME})
target_link_libraries(${PROJECT_NAME} effolkronium_random)
target_include_directories( ${PROJECT_NAME}
    PUBLIC ${PROJECT_SOURCE_DIR}/include
)