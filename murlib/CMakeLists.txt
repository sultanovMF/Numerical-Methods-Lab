project (murlib)

add_library(${PROJECT_NAME} "src/murlib.cpp" "include/murlib/matrix_solver.h" "include/murlib/approximation.h" "include/murlib/integration.h")
add_library(murlib::interpolation ALIAS ${PROJECT_NAME})

target_include_directories( ${PROJECT_NAME}
    PUBLIC ${PROJECT_SOURCE_DIR}/include
)